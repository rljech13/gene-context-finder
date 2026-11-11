"""Read explicit FASTA provenance without guessing source genes from similarity."""
import re
from .sequence import open_text, sha256

ACCESSION=re.compile(r'(?<![A-Za-z0-9])(?:[NXY][PMR]_\d+\.\d+)(?![A-Za-z0-9.])')


def fasta_metadata(path):
    result={}
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith('>'):continue
            header=line[1:].strip();words=header.split()
            if not words:raise ValueError('Empty FASTA header')
            key=words[0]
            if key in result:raise ValueError('Duplicate FASTA identifier: '+key)
            fields=dict(source_header=header);brackets=re.findall(r'\[([^\]]+)\]',header)
            organisms=[x.split('=',1)[1].strip() for x in brackets if x.lower().startswith('organism=')]
            if not organisms:organisms=[x for x in brackets if re.fullmatch(r'[A-Z][a-z]+ [a-z][a-z0-9.-]+(?: [a-z][a-z0-9.-]+)*',x)]
            if len(set(organisms))==1:fields['species']=organisms[0]
            for field,pattern in [('gene_id',r'\bGeneID[=:](\d+)\b'),('gene_symbol',r'\[gene=([^\]]+)\]')]:
                values=re.findall(pattern,header)
                if len(set(values))==1:fields[field]=values[0]
            accessions=set(ACCESSION.findall(header))
            for field,pattern in [('protein',r'^[NXY]P_'),('transcript',r'^[NXY][MR]_')]:
                values=[a for a in accessions if re.match(pattern,a)]
                if len(values)==1:fields[field]=values[0]
            result[key]=fields
    return result


def restore_headers(data,path):
    """Explicit migration of an old checkpoint, only from its unchanged FASTA."""
    digest=sha256(path)
    if digest!=data['query_provenance']['sha256']:raise ValueError('Source-query FASTA differs from the original search input')
    headers=fasta_metadata(path)
    members=[m for q in data['queries'] for m in q['members']]
    if set(headers)!={m['query'] for m in members}:raise ValueError('Source-query identifiers differ from the saved search')
    for q in data['queries']:
        q['members']=[dict(headers[m['query']],**m) for m in q['members']]
    data['query_provenance']['header_metadata']=dict(fasta_sha256=digest,parser='explicit_fasta_fields_v1')
    qmap={q['id']:q for q in data['queries']}
    for locus in data['loci']:
        species=sorted({m['species'] for qid in locus['query_ids'] for m in qmap[qid]['members'] if m.get('species')})
        locus.update(species=species,support_species=len(species))
