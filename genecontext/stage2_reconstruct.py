"""Project a supplied full transcript onto a reviewed genomic window.

Reference-derived RNA is not independent experimental validation. An existing
annotation conflict is resolved only when the reference record also reports
RNA-seq support for all introns. Scientific model sequences are never replaced.
"""
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET

from .sequence import read_fasta, write_fasta, sha256, reverse_complement, translate


def transcripts(path):
    records = {}
    for entry in ET.parse(path).getroot().iter('GBSeq'):
        accession = entry.findtext('GBSeq_accession-version')
        sequence = (entry.findtext('GBSeq_sequence') or '').upper()
        features = entry.findall('GBSeq_feature-table/GBFeature')
        coding = [f for f in features if f.findtext('GBFeature_key') == 'CDS']
        if not accession or not sequence or len(coding) != 1:
            raise ValueError('Each reference record needs an accession, sequence and one CDS')
        location = coding[0].findtext('GBFeature_location') or ''
        match = re.fullmatch(r'(\d+)\.\.(\d+)', location)
        if not match:
            raise ValueError('Only complete, simple transcript CDS intervals are supported: ' + accession)
        start, end = int(match[1])-1, int(match[2])
        qualifiers = {q.findtext('GBQualifier_name'): q.findtext('GBQualifier_value')
                      for q in coding[0].findall('GBFeature_quals/GBQualifier')}
        if qualifiers.get('codon_start', '1') != '1' or not 0 <= start < end <= len(sequence):
            raise ValueError('Invalid/partial transcript CDS: ' + accession)
        notes = [q.findtext('GBQualifier_value') or '' for f in features
                 if f.findtext('GBFeature_key') == 'gene'
                 for q in f.findall('GBFeature_quals/GBQualifier')
                 if q.findtext('GBQualifier_name') == 'note']
        record = dict(accession=accession, sequence=sequence, cds_start=start, cds_end=end,
            protein_accession=qualifiers.get('protein_id', ''), definition=entry.findtext('GBSeq_definition'),
            reference_derived='predicted' in (entry.findtext('GBSeq_definition') or '').lower(),
            source_notes=notes, rna_all_introns_reported=any('rnaseq' in n.lower() and
                'support for all annotated introns' in n.lower() for n in notes))
        if accession in records:
            raise ValueError('Duplicate transcript accession: ' + accession)
        records[accession] = record
    if not records:
        raise ValueError('No GenBank XML transcripts supplied')
    return records


def project_sam(text, transcript, reference_id, dna, native_offset):
    """Require unique, exact full-RNA mapping, and project its CDS, on either strand."""
    records = []
    for line in text.splitlines():
        if not line or line.startswith('@'):
            continue
        f = line.split('\t')
        if len(f) < 11 or f[0] != transcript['accession']:
            continue
        flag = int(f[1])
        if flag & (4 | 256 | 2048):
            continue
        records.append(f)
    if len(records) != 1:
        raise ValueError('Expected one primary transcript alignment')
    f = records[0]; flag = int(f[1]); query = transcript['sequence']
    if f[2] != reference_id or int(f[4]) < 20:
        raise ValueError('Transcript alignment is ambiguous or names another window')
    cigar = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
    if ''.join(n+op for n, op in cigar) != f[5] or not cigar:
        raise ValueError('Invalid transcript CIGAR')
    tags = {p.split(':', 2)[0]: p.split(':', 2)[-1] for p in f[11:] if ':' in p}
    if tags.get('NM') != '0' or any(op not in {'M', '=', 'N'} for _, op in cigar):
        raise ValueError('Full transcript must map exactly, without clips, indels or mismatches')
    reverse = bool(flag & 16); strand = '-' if reverse else '+'
    target = int(f[3])-1; qpos = 0; exons = []; mapped = []; blocks = []
    lo, hi = transcript['cds_start'], transcript['cds_end']
    if reverse:
        lo, hi = len(query)-hi, len(query)-lo
    for number, op in cigar:
        count = int(number)
        if op == 'N':
            target += count
            continue
        if target < 0 or target+count > len(dna):
            raise ValueError('Transcript maps outside the saved DNA window')
        mapped.append(dna[target:target+count].upper()); exons.append([target+native_offset, target+count+native_offset])
        left, right = max(qpos, lo), min(qpos+count, hi)
        if left < right:
            blocks.append([target+left-qpos+native_offset, target+right-qpos+native_offset])
        qpos += count; target += count
    aligned = ''.join(mapped)
    expected = reverse_complement(query) if reverse else query
    if qpos != len(query) or aligned != expected:
        raise ValueError('Reference RNA and mapped genomic exon sequence differ')
    blocks.sort(reverse=reverse); used = 0; phased = []
    for a, b in blocks:
        phased.append([a, b, (3-used % 3) % 3]); used += b-a
    cds = query[transcript['cds_start']:transcript['cds_end']]
    protein = translate(cds)
    if not cds.startswith('ATG') or len(cds) % 3 or not protein.endswith('*') or '*' in protein[:-1]:
        raise ValueError('Reconstructed transcript CDS is not an intact ORF')
    return dict(strand=strand, cds_blocks=phased, exon_blocks=exons,
                cds=cds, protein=protein[:-1], cigar=f[5], mapping_quality=int(f[4]),
                transcript_length=len(query), matched_bases=len(query), sequence_identity=1.0)


def load(root, data, evidence_name='evidence.json', supplement_name='reconstructions.json'):
    path = Path(root) / supplement_name
    if not path.exists():
        return data
    saved = json.loads(path.read_text())
    if (saved.get('schema') != 'stage2_reconstructions_v1' or
            saved.get('source_evidence_sha256') != sha256(Path(root) / evidence_name)):
        raise ValueError('Transcript reconstructions belong to another checkpoint')
    models = {m['model_id']: m for m in data['models']}
    result = copy.deepcopy(data)
    attached = {}
    for rec in saved['models']:
        m = models.get(rec['model_id'])
        if not m or m['source'] != 'annotation':
            raise ValueError('Reconstruction must match an existing annotation model')
        if (rec['cds_sha256'] != m['cds_sha256'] or rec['protein_sha256'] != m['protein_sha256'] or
                rec['cds_blocks'] != m['cds_blocks'] or rec['strand'] != m['strand']):
            raise ValueError('Transcript reconstruction differs from the existing model')
        if (rec['matched_bases'] != rec['transcript_length'] or rec['sequence_identity'] != 1.0
                or rec['mapping_quality'] < 20):
            raise ValueError('Incomplete or ambiguous transcript reconstruction')
        if (rec.get('status') != 'existing_full_model_reconstructed' or
                type(rec.get('rna_all_introns_reported')) is not bool or
                rec.get('resolves_boundary_conflict') is not rec['rna_all_introns_reported']):
            raise ValueError('Invalid transcript support classification')
        if rec['model_id'] in attached:
            raise ValueError('Repeated reconstruction for one model')
        attached[rec['model_id']] = rec
    for m in result['models']:
        if m['model_id'] in attached:
            m['reconstruction'] = attached[m['model_id']]
    return result


def run(args):
    from .stage2_choices import verify_checkpoint
    from .stage2 import TimedRunner
    root, data = verify_checkpoint(args.results)
    out = Path(args.out).resolve()
    if out == root or root in out.parents or out in root.parents:
        raise ValueError('Use a separate reconstruction output directory')
    out.mkdir(parents=True, exist_ok=False)
    references = transcripts(args.transcripts)
    with open(args.mapping, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or any(not r.get('model_id') or not r.get('accession') for r in rows):
        raise ValueError('Mapping TSV needs model_id and accession columns')
    if len({r['model_id'] for r in rows}) != len(rows):
        raise ValueError('Choose one reconstruction transcript per model')
    lookup = {m['model_id']: m for m in data['models']}; regions = {r['region_id']: r for r in data['regions']}
    windows = read_fasta(root/'windows.fna'); coding = read_fasta(root/'cds.fna'); proteins = read_fasta(root/'proteins.faa')
    binary = shutil.which(args.minimap2)
    if not binary:
        raise ValueError('minimap2 is required for splice-aware transcript reconstruction')
    probe = subprocess.run([binary, '--version'], capture_output=True, text=True)
    if probe.returncode:
        raise ValueError('minimap2 version probe failed: ' + probe.stderr[-500:])
    version = probe.stdout.strip()
    runner = TimedRunner(out); results = []; fasta_cds = []; fasta_proteins = []; gff = ['##gff-version 3']
    started = time.perf_counter()
    for i, row in enumerate(rows):
        mid = row['model_id']; model = lookup.get(mid); ref = references.get(row['accession'])
        if not model or model['source'] != 'annotation' or not ref:
            raise ValueError('Unknown annotation model or transcript accession')
        region = regions[model['region_id']]; folder = out/f'projection_{i+1:03d}'; folder.mkdir()
        write_fasta(folder/'window.fna', [(region['region_id'], windows[region['region_id']])])
        write_fasta(folder/'transcript.fna', [(ref['accession'], ref['sequence'])])
        runner.run([binary, '-ax', 'splice', '-uf', '-k14', '-G', str(args.max_intron),
            '-t', args.threads, folder/'window.fna', folder/'transcript.fna'], folder/'alignment.sam')
        projected = project_sam((folder/'alignment.sam').read_text(), ref, region['region_id'],
                                windows[region['region_id']], region['window_start'])
        if (projected['strand'] != model['strand'] or projected['cds_blocks'] != model['cds_blocks'] or
                projected['cds'] != coding[mid] or projected['protein'] != proteins[mid]):
            raise ValueError('Projection proposes a different CDS; retain it for a new-model review, not an automatic replacement')
        rec = {k:v for k,v in projected.items() if k not in {'cds','protein'}}
        rec.update(model_id=mid, accession=ref['accession'], protein_accession=ref['protein_accession'],
            cds_sha256=model['cds_sha256'], protein_sha256=model['protein_sha256'],
            cds_length=len(projected['cds']), protein_length=len(projected['protein']),
            source_notes=ref['source_notes'], reference_derived=ref['reference_derived'],
            rna_all_introns_reported=ref['rna_all_introns_reported'],
            resolves_boundary_conflict=ref['rna_all_introns_reported'],
            status='existing_full_model_reconstructed', alignment_sha256=sha256(folder/'alignment.sam'),
            note='Exact splice projection confirms model consistency. RNA-seq support is reported by the supplied RefSeq record; raw RNA-seq was not reanalysed.')
        header=mid+' transcript='+ref['accession']+' genome='+model['genome_id']+' target_family='+region.get('gene','target')
        if ref['protein_accession']:
            header+=' reference_protein='+ref['protein_accession']
        results.append(rec); fasta_cds.append((header, projected['cds']))
        fasta_proteins.append((header, projected['protein']))
        start=min(a for a,b in rec['exon_blocks']); end=max(b for a,b in rec['exon_blocks'])
        base=f"{model['contig']}\ttranscript_projection\t"
        gff.extend([f"{base}gene\t{start+1}\t{end}\t.\t{model['strand']}\t.\tID={mid}.reconstructed",
                    f"{base}mRNA\t{start+1}\t{end}\t.\t{model['strand']}\t.\tID={mid};Parent={mid}.reconstructed;transcript={ref['accession']}"])
        for a,b in rec['exon_blocks']:
            gff.append(f"{base}exon\t{a+1}\t{b}\t.\t{model['strand']}\t.\tParent={mid}")
        for a,b,p in rec['cds_blocks']:
            gff.append(f"{base}CDS\t{a+1}\t{b}\t.\t{model['strand']}\t{p}\tParent={mid};protein_id={mid}")
    write_fasta(out/'reconstructed.cds.fna', fasta_cds);write_fasta(out/'reconstructed.proteins.faa', fasta_proteins)
    (out/'reconstructed.gff3').write_text('\n'.join(gff)+'\n')
    payload=dict(schema='stage2_reconstructions_v1',source_evidence_sha256=sha256(root/'evidence.json'),
        models=results,provenance=dict(transcripts_sha256=sha256(args.transcripts),mapping_sha256=sha256(args.mapping),
        minimap2_version=version,minimap2_sha256=sha256(binary),commands=runner.commands,
        elapsed_seconds=time.perf_counter()-started,implementation_sha256=sha256(Path(__file__))))
    (out/'reconstructions.json').write_text(json.dumps(payload,indent=2)+'\n')
    if args.attach:
        if (root/'reconstructions.json').exists():
            previous=json.loads((root/'reconstructions.json').read_text())
            signature=lambda record: (record['model_id'],record['accession'],record['cds_sha256'],record['protein_sha256'],record['cds_blocks'])
            if (sorted(map(signature,previous['models']))!=sorted(map(signature,results)) or
                    previous['provenance']['transcripts_sha256']!=payload['provenance']['transcripts_sha256'] or
                    previous['provenance']['mapping_sha256']!=payload['provenance']['mapping_sha256']):
                raise ValueError('Different reconstruction support is already attached; use a separate reviewed revision')
        seal=json.loads((root/'seal.json').read_text());(out/'parent.seal.json').write_text(json.dumps(seal,indent=2)+'\n')
        target=root/'reconstruction';target.mkdir(exist_ok=(root/'reconstructions.json').exists())
        for name in ['reconstructed.cds.fna','reconstructed.proteins.faa','reconstructed.gff3']:
            shutil.copyfile(out/name,target/name);seal['reconstruction/'+name]=sha256(target/name)
        shutil.copyfile(out/'reconstructions.json',root/'reconstructions.json')
        enriched=load(root,data);seal['reconstructions.json']=sha256(root/'reconstructions.json')
        from .stage2_report import write_report
        write_report(root,enriched)
        for name in ['report.html','gene_identities.tsv']:
            seal[name]=sha256(root/name)
        pending=root/'seal.reconstruction.tmp';pending.write_text(json.dumps(seal,indent=2)+'\n');os.replace(pending,root/'seal.json')
    print(f'Reconstructed {len(results)} full model(s): {out}. No review decisions were written.')
    return 0
