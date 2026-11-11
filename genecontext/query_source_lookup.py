"""Verified RefSeq provenance, independent of gene names and genome-search scores."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from .early_synteny import cache_lock

PROTEIN=re.compile(r'(?<![A-Za-z0-9])(?:[NXY]P_\d+\.\d+)(?![A-Za-z0-9.])')
BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'


def record(accession,cache,mode):
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True)
    path=cache/(accession+'.xml');receipt=cache/(accession+'.json')
    url=BASE+'?'+urlencode(dict(db='protein',id=accession,rettype='gb',retmode='xml',tool='gene_context_finder'))
    with cache_lock(cache/'request.lock'):
        if path.exists() or receipt.exists():
            if not path.exists() or not receipt.exists():raise ValueError('Incomplete source metadata cache: '+accession)
            raw=path.read_bytes();info=json.loads(receipt.read_text())
            if info.get('accession')!=accession or hashlib.sha256(raw).hexdigest()!=info.get('sha256'):raise ValueError('Source metadata cache checksum differs: '+accession)
        else:
            if mode=='off':return None
            for attempt in range(2):
                try:
                    with urlopen(Request(url,headers={'User-Agent':'gene-context-finder/source-metadata'}),timeout=20) as response:raw=response.read(4*1024*1024+1)
                    if len(raw)>4*1024*1024:raise ValueError('Source protein metadata exceeds 4 MiB')
                    break
                except OSError:
                    if attempt:raise
                    time.sleep(1)
            # A versioned request must return that exact record before caching.
            root=ET.fromstring(raw)
            if len([r for r in root.iter('GBSeq') if r.findtext('GBSeq_accession-version')==accession])!=1:raise ValueError('NCBI did not return the exact requested accession.version: '+accession)
            info=dict(accession=accession,url=url,sha256=hashlib.sha256(raw).hexdigest(),retrieved_utc=datetime.now(timezone.utc).isoformat())
            path.with_suffix('.xml.part').write_bytes(raw);path.with_suffix('.xml.part').replace(path)
            receipt.write_text(json.dumps(info,indent=2)+'\n');time.sleep(.35)
    root=ET.fromstring(raw);items=[r for r in root.iter('GBSeq') if r.findtext('GBSeq_accession-version')==accession]
    if len(items)!=1:raise ValueError('Cached accession.version differs: '+accession)
    item=items[0];values={}
    for q in item.iter('GBQualifier'):values.setdefault(q.findtext('GBQualifier_name'),set()).add(q.findtext('GBQualifier_value') or '')
    genes=sorted({v.split(':',1)[1] for v in values.get('db_xref',[]) if re.fullmatch(r'GeneID:\d+',v)})
    symbols=sorted(values.get('gene',[]))
    return dict(info,sequence=''.join((item.findtext('GBSeq_sequence') or '').split()).upper(),species=item.findtext('GBSeq_organism') or '',gene_ids=genes,gene_symbol=symbols[0] if len(symbols)==1 else '')


def enrich(queries,cache,out,mode='auto'):
    """Only verified sequence + one GeneID may supply otherwise missing metadata."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=[];records={}
    for query in queries:
        for member in query['members']:
            if member.get('gene_id'):continue
            accessions={a for field in ['query','protein'] for a in PROTEIN.findall(member.get(field,''))}
            if len(accessions)!=1:continue
            accession=next(iter(accessions));row=dict(query=member['query'],accession=accession,status='unresolved')
            try:
                if accession not in records:records[accession]=record(accession,cache,mode)
                ref=records[accession]
                if ref:
                    target=out/'records'/(accession+'.xml');target.parent.mkdir(exist_ok=True)
                    shutil.copyfile(Path(cache)/(accession+'.xml'),target)
                    if hashlib.sha256(target.read_bytes()).hexdigest()!=ref['sha256']:raise ValueError('Source metadata changed while copying')
                    row['record_file']='records/'+target.name
                if ref is None:row['status']='lookup_disabled_cache_missing'
                elif ref['sequence']!=query['sequence']:row.update(status='query_differs_from_accession_sequence',source_sequence_sha256=hashlib.sha256(ref['sequence'].encode()).hexdigest())
                elif member.get('species') and member['species']!=ref['species']:row['status']='source_species_conflict'
                elif len(ref['gene_ids'])!=1:row['status']='source_gene_id_ambiguous_or_missing'
                else:
                    member['gene_id']=ref['gene_ids'][0]
                    for field in ['species','gene_symbol']:
                        if not member.get(field) and ref[field]:member[field]=ref[field]
                    member['source_metadata']=dict(accession=accession,url=ref['url'],xml_sha256=ref['sha256'],sequence_verified=True,basis='exact_versioned_refseq_record')
                    row.update(status='verified_refseq_metadata',gene_id=member['gene_id'])
                if ref:row.update(url=ref['url'],xml_sha256=ref['sha256'])
            except (OSError,ValueError,RuntimeError,ET.ParseError) as error:row.update(status='metadata_lookup_unavailable',error=str(error))
            rows.append(row)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'receipt.json').write_text(json.dumps(dict(mode=mode,records=rows,sequence_replacement=False),indent=2)+'\n')
    return dict(mode=mode,records=rows,sequence_replacement=False)
