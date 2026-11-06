#!/usr/bin/env python3
"""Fetch a versioned RefSeq POU5 test panel (stdlib only; no inferred sequences).

This is a named-gene test panel, not an assertion that POU5F1 and all POU5F3
copies are one-to-one orthologs. All NM_/XM_ products listed by these NCBI Gene
records are considered, across their reported assemblies. Raw XML is retained.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import gzip
import json
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genecontext.sequence import reverse_complement, translate, write_fasta, write_tsv

GENES = ['5460', '30333', '427781', '100498076', '100170181', '100170593']


def fetch(root, name, **params):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(params)
    path = root/'sources'/(name+'.gz')
    if not path.exists():
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    content = response.read()
                if not content or b'<ERROR>' in content:
                    raise ValueError('NCBI returned an empty/error response')
                path.write_bytes(gzip.compress(content,mtime=0))
                break
            except (OSError, ValueError):
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
        time.sleep(0.4)
    return path, url


def write_table(path, records, columns):
    write_tsv(path, records, columns)
    path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))


def qualifiers(feature):
    return {q.findtext('GBQualifier_name'): q.findtext('GBQualifier_value', '')
            for q in feature.findall('./GBFeature_quals/GBQualifier')}


def build(output, controls=False):
    root = Path(output).resolve(); (root/'sources').mkdir(parents=True, exist_ok=True)
    source, url = fetch(root, 'genes.xml', db='gene', id=','.join(GENES), retmode='xml')
    provenance = [dict(file=str(source.relative_to(root)), url=url)]
    records, primary = {}, []
    for gene in ET.fromstring(gzip.decompress(source.read_bytes())).findall('Entrezgene'):
        gid = gene.findtext('.//Gene-track_geneid')
        symbol, species = gene.findtext('.//Gene-ref_locus'), gene.findtext('.//Org-ref_taxname')
        genomic = gene.findall('./Entrezgene_locus/Gene-commentary')
        preferred = next((g for g in genomic if 'Primary Assembly' in g.findtext('Gene-commentary_heading', '')), genomic[0])
        seq = preferred.find('./Gene-commentary_seqs//Seq-interval')
        if seq is not None:
            primary.append(dict(gene_id=gid, symbol=symbol, species=species,
                                contig=preferred.findtext('Gene-commentary_accession')+'.'+preferred.findtext('Gene-commentary_version'),
                                start=int(seq.findtext('Seq-interval_from')), end=int(seq.findtext('Seq-interval_to'))+1,
                                strand='-' if seq.find('Seq-interval_strand/Na-strand').get('value') == 'minus' else '+',
                                assembly=preferred.findtext('Gene-commentary_heading')))
        for product in gene.findall('.//Gene-commentary'):
            accession = product.findtext('Gene-commentary_accession', '')
            if accession.startswith(('NM_', 'XM_')):
                version = accession+'.'+product.findtext('Gene-commentary_version')
                records[version] = dict(gene_id=gid, gene_symbol=symbol, species=species)
    proteins, cds_records, rows, excluded = [], [], [], []
    accessions = sorted(records)
    for begin in range(0, len(accessions), 15):
        source, url = fetch(root, f'transcripts_{begin//15:02d}.xml', db='nuccore',
                            id=','.join(accessions[begin:begin+15]), rettype='gb', retmode='xml')
        provenance.append(dict(file=str(source.relative_to(root)), url=url))
        for record in ET.fromstring(gzip.decompress(source.read_bytes())).findall('GBSeq'):
            accession = record.findtext('GBSeq_accession-version')
            if accession not in records:
                raise ValueError('NCBI returned a different accession/version: '+str(accession))
            seq = record.findtext('GBSeq_sequence').upper()
            features = [f for f in record.findall('./GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key') == 'CDS']
            for feature in features:
                q = qualifiers(feature)
                if 'translation' not in q or 'pseudo' in q:
                    excluded.append(dict(transcript=accession, reason='no translated CDS')); continue
                if 'GeneID:'+records[accession]['gene_id'] not in [x.findtext('GBQualifier_value') for x in feature.findall('./GBFeature_quals/GBQualifier')]:
                    raise ValueError('CDS GeneID mismatch: '+accession)
                pieces = []
                for interval in feature.findall('./GBFeature_intervals/GBInterval'):
                    if interval.findtext('GBInterval_accession', accession) != accession:
                        raise ValueError('Remote CDS interval: '+accession)
                    a, b = int(interval.findtext('GBInterval_from'))-1, int(interval.findtext('GBInterval_to'))
                    part = seq[min(a,b-1):max(a+1,b)]
                    if interval.find('GBInterval_iscomp') is not None and interval.find('GBInterval_iscomp').get('value') == 'true':
                        part = reverse_complement(part)
                    pieces.append(part)
                cds = ''.join(pieces)[int(q.get('codon_start', '1'))-1:]
                protein = ''.join(q['translation'].split())
                translated = translate(cds).removesuffix('*')
                start_exception = bool(translated and translated[1:] == protein[1:] and protein.startswith('M') and translated != protein)
                if translated != protein and not start_exception:
                    raise ValueError('CDS/annotated translation mismatch: '+accession)
                if len(cds)%3 or '*' in translated or len(protein)<30:
                    excluded.append(dict(transcript=accession, reason='invalid/short translated CDS')); continue
                metadata = records[accession]
                qid = metadata['species'].replace(' ', '_')+'__'+metadata['gene_symbol']+'__'+accession
                proteins.append((qid, protein)); cds_records.append((qid, cds))
                rows.append(dict(query=qid, group='POU5', **metadata, transcript=accession,
                                 protein=q['protein_id'], annotation_status='predicted RefSeq' if accession.startswith('XM_') else 'non-predicted RefSeq',
                                 protein_aa=len(protein), cds_bp=len(cds), initiation_exception=start_exception,
                                 gene_url='https://www.ncbi.nlm.nih.gov/gene/'+metadata['gene_id'],
                                 transcript_url='https://www.ncbi.nlm.nih.gov/nuccore/'+accession,
                                 protein_url='https://www.ncbi.nlm.nih.gov/protein/'+q['protein_id'],
                                 protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
                                 cds_sha256=hashlib.sha256(cds.encode()).hexdigest()))
    missing = set(accessions)-{r['transcript'] for r in rows}-{r['transcript'] for r in excluded}
    if missing:
        raise ValueError('Missing fetched transcripts: '+str(sorted(missing)))
    write_fasta(root/'proteins.faa', proteins); write_fasta(root/'cds.fna', cds_records)
    write_table(root/'queries.tsv', rows, list(rows[0]))
    write_table(root/'reference_loci.tsv', primary, list(primary[0]))
    if controls:
        genomes = []; (root/'controls').mkdir(exist_ok=True)
        for item in primary:
            start, end = max(0, item['start']-10000), item['end']+10000
            name = item['species'].replace(' ', '_')+'__'+item['symbol']
            source, url = fetch(root, name+'.fna', db='nuccore', id=item['contig'],
                                seq_start=start+1, seq_stop=end, rettype='fasta', retmode='text')
            provenance.append(dict(file=str(source.relative_to(root)), url=url))
            # Keep original genomic accession in the source header; use a clear
            # local window ID in the control FASTA and record the offset.
            lines = gzip.decompress(source.read_bytes()).decode().splitlines(); sequence = ''.join(lines[1:])
            if not lines[0].startswith('>') or len(sequence) != end-start:
                raise ValueError('Unexpected control window: '+name)
            target = root/'controls'/(name+'.fna'); write_fasta(target, [(name, sequence)])
            genomes.append(dict(genome_id=name, fasta=str(target.relative_to(root)), gff3='',
                                species=item['species'], original_contig=item['contig'], source_start=start,
                                expected_start=item['start']-start, expected_end=item['end']-start,
                                expected_strand=item['strand']))
        write_table(root/'controls.genomes.tsv', genomes, list(genomes[0]))
    for entry in provenance:
        entry['sha256'] = hashlib.sha256((root/entry['file']).read_bytes()).hexdigest()
        entry['retrieved_sha256'] = hashlib.sha256(gzip.decompress((root/entry['file']).read_bytes())).hexdigest()
    (root/'provenance.json').write_text(json.dumps(dict(created=datetime.now(timezone.utc).isoformat(),
        source='NCBI Gene and RefSeq annotated transcript CDS', genes=GENES,
        scope='All NM_/XM_ products listed in the fetched Gene records; not all annotation providers or releases',
        caution='POU5F1 and POU5F3 copies are labelled separately; group POU5 is a search panel, not an orthology assignment',
        queries=len(rows), excluded=excluded, sources=provenance), indent=2)+'\n')
    print(json.dumps(dict(output=str(root), variants=len(rows), species=sorted({r['species'] for r in rows}), excluded=excluded), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    parser.add_argument('--controls', action='store_true', help='Also download small genomic reference windows for smoke tests')
    args = parser.parse_args(); build(args.out, args.controls)
