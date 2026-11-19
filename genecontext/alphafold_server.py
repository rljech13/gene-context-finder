"""Read-only AlphaFold Server result adapter; requests and old results stay sealed.

Provider-specific ZIP/mmCIF details end here. STOP4 consumes the ordinary exact
sequence catalogue, so future providers can produce that same boundary object.
"""
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
import zipfile

from .sequence import sha256
from .stage1 import dump, new_output, now
from .stage4_structure import ca_atoms, validate_pae, RESIDUES


def request_protein(job):
    if isinstance(job,list):
        if len(job)!=1:raise ValueError('Expected one returned Server job')
        job=job[0]
    if not isinstance(job,dict) or job.get('dialect')!='alphafoldserver':
        raise ValueError('Expected an AlphaFold Server job request')
    entities=job.get('sequences',[])
    if len(entities)!=1 or set(entities[0])!={'proteinChain'}:
        raise ValueError('Server prototype imports one unmodified protein monomer')
    protein=entities[0]['proteinChain']
    if protein.get('count')!=1 or protein.get('glycans') or protein.get('modifications'):
        raise ValueError('Server prototype imports one unmodified protein monomer')
    sequence=protein.get('sequence','')
    if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('Noncanonical returned Server protein')
    return job,sequence


def convert_monomer(cif_text, output, sequence):
    try:import gemmi
    except ImportError as error:
        raise ValueError('AF3 mmCIF import needs the optional structures extra: pip install "gene-context-finder[structures]"') from error
    structure=gemmi.make_structure_from_block(gemmi.cif.read_string(cif_text).sole_block())
    if len(structure)!=1 or len(structure[0])!=1:
        raise ValueError('Expected a single model and single protein chain')
    chain=structure[0][0];original_chain=chain.name
    residues=[]
    for residue in chain:
        if residue.name not in RESIDUES:raise ValueError('Noncanonical mmCIF residue')
        residues.append(RESIDUES[residue.name])
        for atom in residue:
            if not all(math.isfinite(v) for v in [atom.pos.x,atom.pos.y,atom.pos.z,atom.b_iso]) or not 0<=atom.b_iso<=100:
                raise ValueError('Invalid AF3 coordinate or pLDDT')
    if ''.join(residues)!=sequence:raise ValueError('mmCIF is not the exact requested full protein')
    chain.name='A'
    output.write_text(structure.make_pdb_string())
    atoms=ca_atoms(output,'A')
    if ''.join(a['aa'] for a in atoms)!=sequence:raise ValueError('Converted PDB lacks exact full CA coverage')
    return dict(converter='gemmi',converter_version=gemmi.__version__,original_chain=original_chain,
        output_chain='A',note='Raw mmCIF preserved in the original ZIP; PDB coordinates/confidence have PDB decimal precision.')


def import_result(requests, archive, out, source_url=''):
    request_file=Path(requests).resolve();plan=json.loads(request_file.read_text())
    if plan.get('schema')!='gcf_af3_requests_v1':raise ValueError('Expected prepared gcf_af3_requests_v1 requests')
    # Validate inputs before creating output. Retain source files verbatim.
    for request in plan['requests']:
        for key in ['input','server_input']:
            if sha256(request_file.parent/request[key+'_file'])!=request[key+'_sha256']:
                raise ValueError('Prepared AF3 request changed')
    archive=Path(archive).resolve()
    if archive.stat().st_size>256*1024**2:raise ValueError('AF3 archive exceeds prototype size limit')
    with zipfile.ZipFile(archive) as zipped:
        names=zipped.namelist()
        if len(names)!=len(set(names)) or len(names)>2000 or sum(f.file_size for f in zipped.infolist())>512*1024**2:
            raise ValueError('Unsupported AF3 archive size or duplicate entries')
        if any(PurePosixPath(n).is_absolute() or '..' in PurePosixPath(n).parts or '\\' in n for n in names):
            raise ValueError('Unsafe AF3 archive member')
        def read(name):
            if zipped.getinfo(name).file_size>64*1024**2:raise ValueError('AF3 result member exceeds size limit')
            return zipped.read(name)
        job_files=[n for n in names if n.endswith('_job_request.json')]
        if len(job_files)!=1:raise ValueError('Expected one Server job request in archive')
        job,sequence=request_protein(json.loads(read(job_files[0])))
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        matches=[r for r in plan['requests'] if r['protein_sha256']==digest]
        if len(matches)!=1:raise ValueError('Returned sequence does not match a prepared protein')
        request=matches[0]
        expected,expected_sequence=request_protein(json.loads((request_file.parent/request['server_input_file']).read_text()))
        if sequence!=expected_sequence:raise ValueError('Returned protein differs from original prepared input')
        if job['name']!=expected['name']:raise ValueError('Returned job name differs from prepared request')
        if list(map(str,job.get('modelSeeds',[])))!=list(map(str,expected.get('modelSeeds',[]))):
            raise ValueError('Returned seed differs from prepared request')
        candidates=[]
        for name in names:
            match=re.fullmatch(r'(.+)_summary_confidences(_[0-9]+)?\.json',name)
            if not match:continue
            stem,index=match[1],match[2] or ''
            model=stem+'_model'+index+'.cif'
            if model not in names:continue
            summary=json.loads(read(name));score=summary.get('ranking_score')
            if type(score) not in {int,float} or not math.isfinite(score):raise ValueError('Invalid AF3 ranking score')
            candidates.append((score,model,stem+'_full_data'+index+'.json',summary))
        if not candidates:raise ValueError('No ranked Server mmCIF model found')
        score,model,full,summary=max(candidates,key=lambda r:(r[0],r[1]))
        root=new_output(out)
        try:
            original=root/'source.zip';shutil.copyfile(archive,original)
            cif=root/'model.cif';cif.write_bytes(read(model))
            pdb=root/'model.pdb';conversion=convert_monomer(cif.read_text(),pdb,sequence)
            dump(root/'summary_confidences.json',summary)
            (root/'job_request.json').write_bytes(read(job_files[0]))
            if 'terms_of_use.md' in names:(root/'terms_of_use.md').write_bytes(read('terms_of_use.md'))
            (root/'README.txt').write_text('AlphaFold Server output: https://alphafoldserver.com/output-terms\nNon-commercial use only. No use in docking/screening tools or training biomolecular structure prediction models.\nAttribution: Abramson et al., Nature (2024), doi:10.1038/s41586-024-07487-w.\nModifications: highest-ranked sample selected; original mmCIF converted to chain-A PDB with Gemmi. Raw results and all other samples remain in source.zip.\n')
            entry=dict(file='model.pdb',file_sha256=sha256(pdb),protein_sha256=digest,chain='A',
                source='AlphaFold Server / AlphaFold 3',accession=job['name'],source_url=source_url,
                confidence='plddt_bfactor',provider='alphafold_server',model='AlphaFold 3',
                usage_terms_url='https://alphafoldserver.com/output-terms',
                adapter_version=1,source_archive_sha256=sha256(original),raw_mmcif_sha256=sha256(cif))
            pae_status='not_available'
            if full in names:
                raw=read(full);(root/'full_data.json').write_bytes(raw);data=json.loads(raw)
                matrix=data.get('pae',data.get('full_pae'))
                if matrix is not None:
                    dump(root/'pae.json',dict(predicted_aligned_error=matrix));validate_pae(root/'pae.json',len(sequence))
                    entry.update(pae_file='pae.json',pae_file_sha256=sha256(root/'pae.json'));pae_status='available_not_interpreted'
            entries=[dict(entry,sequence_id=sid) for sid in request['source_sequence_ids']]
            dump(root/'catalogue.json',dict(schema='stage4_structures_v1',entries=entries))
            dump(root/'receipt.json',dict(schema='gcf_structure_result_v1',provider='alphafold_server',model='AlphaFold 3',
                adapter_version=1,created=now(),status='available',job_name=job['name'],job_source_url=source_url,
                server_request_version=job.get('version'),model_seeds=job.get('modelSeeds'),
                protein_sha256=digest,source_sequence_ids=request['source_sequence_ids'],
                request_manifest_sha256=sha256(request_file),source_archive_sha256=sha256(original),
                selected_member=model,ranking_score=score,returned_samples=len(candidates),
                selection_policy='Highest returned ranking_score, filename breaks ties; all samples retained in source.zip',
                raw_mmcif_sha256=sha256(cif),pae_status=pae_status,conversion=conversion,extensions={},
                source_code_sha256=sha256(Path(__file__)),limitations=['A predicted structure does not establish gene identity, orthology or biological function.']))
            dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
        except Exception:
            # Failed attempts are visible and never mistaken for ready outputs.
            dump(root/'import.failed.json',dict(status='failed',source_archive_sha256=sha256(archive)))
            raise
    return root/'catalogue.json'
