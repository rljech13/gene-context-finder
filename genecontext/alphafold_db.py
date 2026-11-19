"""Optional AFDB retrieval for exact selected proteins; never submit folding jobs.

The MD5 lookup is the public AFDB sequence endpoint. SHA-256 remains our storage
identity. API metadata and actual PDB-chain sequence must both match in full.
"""
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import time
from urllib.error import HTTPError
from urllib.parse import urlencode, urlsplit, quote
from urllib.request import Request, urlopen
from uuid import uuid4

from .artifact_cache import cached
from .early_synteny import token
from .sequence import sha256, write_tsv
from .stage1 import dump, now, new_output
from .stage4_structure import ca_atoms, validate_pae

BASE = 'https://alphafold.ebi.ac.uk'
DEFAULT_CACHE = str(Path.home()/'.cache/gene-context-finder/afdb')
POLICY = 'afdb_exact_monomer_v1'


def trusted_url(url):
    parsed = urlsplit(url)
    if parsed.scheme != 'https' or parsed.hostname not in {'alphafold.ebi.ac.uk', 'alphafold.com', 'www.alphafold.com'} or parsed.port not in {None, 443} or parsed.username or parsed.password:
        raise ValueError('AFDB returned an unexpected download URL')
    return url


class Download:
    def __init__(self):
        self.last = 0.

    def __call__(self, url, maximum):
        trusted_url(url)
        for attempt in range(3):
            time.sleep(max(0, 1-(time.monotonic()-self.last)))
            self.last = time.monotonic()
            try:
                with urlopen(Request(url, headers={'User-Agent': 'GeneContextFinder/AFDB-prototype', 'Accept': '*/*'}), timeout=25) as response:
                    trusted_url(response.url)
                    raw = response.read(maximum+1)
                    if len(raw)>maximum:
                        raise ValueError('AFDB response exceeds the prototype size limit')
                    return raw
            except HTTPError as error:
                if error.code == 404:
                    return None
                if error.code not in {429, 502, 503, 504} or attempt == 2:
                    raise
                time.sleep(2**(attempt+1))


class Client:
    def __init__(self, cache, *, offline=False, refresh=False, transport=None):
        if offline and refresh:
            raise ValueError('Offline retrieval cannot refresh the AFDB cache')
        self.cache = Path(cache).expanduser().resolve()
        self.offline = offline
        self.transport = transport or Download()
        self.generation = uuid4().hex if refresh else 'pinned'

    def fetch(self, sequence):
        if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Exact structure retrieval requires a canonical full protein sequence')
        if len(sequence)>1500:
            raise ValueError('AFDB prototype limit is 1500 residues; use explicit catalogue import for longer proteins')
        identity = dict(policy=POLICY, protein_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                        implementation_sha256=sha256(Path(__file__)), generation=self.generation)
        expected = self.cache/token(identity)
        if self.offline and not (expected/'cache.json').is_file():
            return None, dict(status='offline_cache_miss', reason='No pinned exact-sequence AFDB response is cached')

        def build(folder):
            request_url = BASE+'/api/sequence/summary?'+urlencode(dict(id=hashlib.md5(sequence.encode()).hexdigest(), type='md5'))
            raw = self.transport(request_url, 4*1024**2)
            result = dict(status='not_found', reason='No full exact-sequence AlphaFold DB monomer found', requested=now(),
                          sequence_sha256=identity['protein_sha256'], lookup_url=request_url, attempts=[])
            if raw is None:
                dump(folder/'result.json', result); return
            (folder/'sequence-summary.json').write_bytes(raw)
            summary = json.loads(raw)
            if not isinstance(summary,dict) or not isinstance(summary.get('entry'),dict) or not isinstance(summary.get('structures'),list):
                raise ValueError('AFDB sequence summary is not an object')
            if summary.get('entry', {}).get('sequence') != sequence:
                raise ValueError('AFDB lookup returned a different full sequence')
            candidates = []
            for item in summary.get('structures', []):
                if not isinstance(item,dict) or not isinstance(item.get('summary'),dict):
                    raise ValueError('AFDB returned an invalid structure summary')
                item = item.get('summary', {})
                mid = item.get('model_identifier', '')
                if item.get('provider') != 'AlphaFold DB' or item.get('coverage') != 1 or item.get('sequence_identity') != 1 or item.get('oligomeric_state') != 'MONOMER':
                    continue
                if not re.fullmatch(r'AF-[A-Za-z0-9-]{1,100}', mid):
                    continue
                score = item.get('confidence_avg_local_score')
                if not isinstance(score, (int,float)) or not math.isfinite(score):
                    score = -1
                candidates.append((item.get('created', ''), score, mid))
            # Deterministic choice among identical full sequences. It is not a
            # species/gene assignment; all source accession metadata is retained.
            for _, _, mid in sorted(candidates, reverse=True)[:3]:
                metadata_url = BASE+'/api/prediction/'+quote(mid, safe='-')
                raw = self.transport(metadata_url, 4*1024**2)
                if raw is None:
                    result['attempts'].append(dict(model=mid, status='metadata_not_found')); continue
                (folder/(mid+'.metadata.json')).write_bytes(raw)
                records = json.loads(raw)
                if not isinstance(records,list) or any(not isinstance(r,dict) for r in records):
                    raise ValueError('AFDB prediction metadata is not a model list')
                accepted = [r for r in records if (r.get('modelEntityId') or r.get('entryId')) == mid
                            and (r.get('sequence') or r.get('uniprotSequence')) == sequence
                            and r.get('providerId') == 'GDM' and not r.get('isComplex', False)
                            and r.get('sequenceStart', r.get('uniprotStart')) == 1
                            and r.get('sequenceEnd', r.get('uniprotEnd')) == len(sequence)]
                if not accepted:
                    result['attempts'].append(dict(model=mid, status='metadata_rejected')); continue
                metadata = accepted[0]; chain = metadata.get('chainId') or 'A'
                if chain != 'A' or not metadata.get('pdbUrl'):
                    result['attempts'].append(dict(model=mid, status='unsupported_coordinate_format')); continue
                raw = self.transport(trusted_url(metadata['pdbUrl']), 20*1024**2)
                if raw is None:
                    result['attempts'].append(dict(model=mid, status='coordinates_not_found')); continue
                pdb = folder/'model.pdb'; pdb.write_bytes(raw)
                atoms = ca_atoms(pdb, chain)
                if ''.join(a['aa'] for a in atoms) != sequence:
                    pdb.unlink(); result['attempts'].append(dict(model=mid, status='coordinate_sequence_mismatch')); continue
                if any(not 0<=a['confidence']<=100 for a in atoms):
                    raise ValueError('AFDB coordinates contain invalid pLDDT')
                result.update(status='available', reason='', model_id=mid, version=metadata.get('latestVersion'),
                    chain=chain, pdb_url=metadata['pdbUrl'], pdb_sha256=sha256(pdb),
                    accession=metadata.get('uniprotAccession') or mid, metadata_url=metadata_url,
                    metadata_file=mid+'.metadata.json', provider='AlphaFold DB', tool=metadata.get('toolUsed'),
                    source_organism=metadata.get('organismScientificName'),
                    selection_policy='Most recent model date, then mean confidence, then ID; only identical full-sequence monomers',
                    mean_plddt=sum(a['confidence'] for a in atoms)/len(atoms), pae_status='not_available')
                if metadata.get('paeDocUrl'):
                    result['pae_url'] = metadata['paeDocUrl']
                    try:
                        raw = self.transport(trusted_url(metadata['paeDocUrl']), 64*1024**2)
                        if raw is not None:
                            pae = folder/'pae.json'; pae.write_bytes(raw); validate_pae(pae, len(sequence))
                            result.update(pae_status='available_not_interpreted', pae_sha256=sha256(pae))
                    except (OSError,ValueError,TypeError,KeyError) as error:
                        (folder/'pae.json').unlink(missing_ok=True)
                        result.update(pae_status='unavailable', pae_error=str(error))
                break
            dump(folder/'result.json', result)

        directory, receipt = cached(self.cache, identity, build)
        result = json.loads((directory/'result.json').read_text())
        return directory, dict(result, cache=receipt)


def collect(sequences, out, cache=DEFAULT_CACHE, limit=3, offline=False, refresh=False, transport=None):
    if not 1<=limit<=50:
        raise ValueError('AFDB prototype limit must be between 1 and 50 unique proteins')
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}', sid) for sid in sequences):
        raise ValueError('Structure sequence IDs must be safe standalone identifiers')
    client = Client(cache, offline=offline, refresh=refresh, transport=transport)
    root = new_output(out); entries=[]; rows=[]; seen={}
    for sid, sequence in sequences.items():
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest not in seen:
            if len(seen)>=limit:
                result = dict(status='not_requested', reason='Pilot unique-protein limit reached'); directory=None
            else:
                try:
                    directory, result = client.fetch(sequence)
                except (OSError,ValueError,TypeError,KeyError) as error:
                    directory=None; result=dict(status='unavailable', reason=str(error))
            seen[digest]=(directory,result)
        directory,result = seen[digest]
        row=dict(sequence_id=sid, protein_sha256=digest, status=result['status'],
                 reason=result.get('reason',''), model_id=result.get('model_id',''),
                 source_organism=result.get('source_organism',''),mean_plddt=result.get('mean_plddt',''),
                 pae_status=result.get('pae_status','not_available'))
        rows.append(row)
        if result['status']!='available':
            continue
        folder=root/'files'/digest;folder.mkdir(parents=True, exist_ok=True)
        for file in directory.iterdir():
            if file.is_file() and file.name!='cache.json' and not (folder/file.name).exists():
                shutil.copyfile(file, folder/file.name)
        entry=dict(sequence_id=sid, file=str((folder/'model.pdb').relative_to(root)),
            file_sha256=result['pdb_sha256'], protein_sha256=digest, chain=result['chain'],
            source='AlphaFold DB', accession=result['accession'], source_url=result['metadata_url'],
            confidence='plddt_bfactor', afdb_model_id=result['model_id'], afdb_version=result['version'],
            source_organism=result.get('source_organism'), pae_status=result['pae_status'])
        if result['pae_status']=='available_not_interpreted':
            entry.update(pae_file=str((folder/'pae.json').relative_to(root)), pae_file_sha256=result['pae_sha256'])
        entries.append(entry)
    dump(root/'catalogue.json',dict(schema='stage4_structures_v1', entries=entries))
    dump(root/'retrieval.json',dict(schema='gcf_afdb_retrieval_v1',provider='AlphaFold DB',policy=POLICY,
        created=now(),rows=rows,requests=[dict(protein_sha256=k,**v[1]) for k,v in seen.items()],
        prediction_submitted=False,source_code_sha256=sha256(Path(__file__)),
        note='Existing predictions only. Identical sequences may share a structure from another accession/organism; this does not assign gene identity or orthology. Missing exact structures require a separate prediction step.'))
    write_tsv(root/'retrieval.tsv', rows, list(rows[0]) if rows else ['sequence_id','status'])
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    return root/'catalogue.json',rows


def add_options(command):
    command.add_argument('--structure-cache',default=DEFAULT_CACHE,help='Shared immutable exact-sequence AFDB download cache')
    command.add_argument('--structure-limit',type=int,default=3,help='Pilot maximum unique proteins to look up (1–50; default 3)')
    command.add_argument('--structure-offline',action='store_true',help='Reuse verified AFDB cache only; no requests')
    command.add_argument('--structure-refresh',action='store_true',help='Explicitly create a new pinned AFDB response; previous cache retained')
