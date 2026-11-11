"""A compact reference sheet derived solely from sealed stage-3 measurements."""
import json
from collections import Counter
from pathlib import Path

from .sequence import write_fasta, write_tsv
from .stage3_analysis import read_table


def oriented_pair(rows, model, reference):
    """Keep both coverage directions when reusing an existing pair table."""
    row = next((p for p in rows if {p['sequence_a'], p['sequence_b']} == {model, reference}), None)
    if row is None:
        return None
    reverse = row['sequence_a'] != model
    return dict(identity=row['identity'], paired_positions=row['paired_positions'],
                identical_positions=row['identical_positions'],
                coverage_model=row['coverage_b' if reverse else 'coverage_a'],
                coverage_reference=row['coverage_a' if reverse else 'coverage_b'])


def sequence_guidance(model):
    quality = model['quality']
    if quality['coding_status'] != 'complete_predicted' or str(model.get('window_edge')).lower() == 'true':
        return dict(label='Check completeness', tone='review',
                    detail='The selected model is partial or touches a saved window edge. Verify its CDS boundaries before use.')
    if quality['database_status'] != 'homology_assignment':
        return dict(label='Check gene identity', tone='review',
                    detail='The coding sequence passes integrity checks, but broad-database naming is unresolved. The original annotation and the selected model are retained.')
    if quality['source_agreement'] == 'different_CDS_models':
        return dict(label='Sources differ', tone='review',
                    detail='Use the selected native sequence as the reviewed candidate. Annotation and Tiberius propose different CDS structures; the chosen accession is retained.')
    return dict(label='Selected native CDS', tone='selected',
                detail='The selected native CDS and protein pass sequence-integrity checks. This is a sequence candidate, not experimental validation or a claim of interchangeability.')


def write_summary(root, payload):
    root = Path(root)
    models = [dict(m) for m in payload['models']]
    counts = Counter((m['group'], m['species']) for m in models)
    numbers = Counter()
    for m in models:
        key = m['group'], m['species']; numbers[key] += 1
        m['copy_number'] = numbers[key] if counts[key] > 1 else None
        m['guidance'] = sequence_guidance(m)
    references, comparisons, sheets = {}, [], []
    for group in payload['groups']:
        gid = group['id']; folder = root/'queries'/gid
        protein = read_table(folder/'protein_pairs.tsv')
        coding = read_table(folder/'cds_pairs.tsv')
        chosen = [m for m in models if m['group'] == gid]
        aligned = payload['alignments'].get(gid, {})
        refs = []
        for q in group['queries']:
            rid = 'reference_'+q
            if rid not in aligned.get('protein', {}):
                continue
            members = payload['query_members'].get(q, []) or [{}]
            first = members[0]
            refs.append(dict(id=rid, kind='protein', query_id=q,
                species=first.get('species') or 'Unspecified species',
                symbol=first.get('gene_symbol') or first.get('reference_label') or group['name'],
                accession=first.get('protein') or first.get('query') or q,
                aliases=members, length=len(aligned['protein'][rid].replace('-', '')),
                origin='Supplied protein query'))
        reference_species = {ref['species'] for ref in refs}
        cds_refs = [dict(id=m['model_id'], kind='cds', species=m['species'], symbol=m['gene_label'],
                        accession=m['original_model'], length=len(aligned['cds'][m['model_id']].replace('-', '')),
                        copy_number=m['copy_number'], origin='Selected STOP-2 model')
                    for m in chosen if m['species'] in reference_species and m['model_id'] in aligned.get('cds', {})]
        references[gid] = dict(protein=refs, cds=cds_refs)
        for m in chosen:
            options = []
            for kind, panel, rows in [('protein', refs, protein), ('cds', cds_refs, coding)]:
                for ref in panel:
                    # Self identity is exact but does not add an independent pair.
                    self_reference = m['model_id'] == ref['id']
                    pair = None if self_reference else oriented_pair(rows, m['model_id'], ref['id'])
                    if not pair and not self_reference:
                        continue
                    row = dict(group=gid, model_id=m['model_id'], species=m['species'],
                               original_model=m['original_model'], reference_id=ref['id'],
                               reference_species=ref['species'], reference_accession=ref['accession'],
                               reference_symbol=ref['symbol'], kind=kind,
                               origin=ref['origin'], self_reference=self_reference,
                               identity='', coverage_model='', coverage_reference='',
                               paired_positions='', identical_positions='', alignment_file='')
                    if self_reference:
                        alphabet = set('ACDEFGHIKLMNPQRSTVWY' if kind == 'protein' else 'ACGT')
                        comparable = sum(c in alphabet for c in aligned[kind][m['model_id']])
                        row.update(paired_positions=comparable, identical_positions=comparable,
                                   identity=1.0 if comparable else '',
                                   coverage_model=1.0 if comparable else '',
                                   coverage_reference=1.0 if comparable else '')
                    if pair:
                        row.update(pair)
                        a, b = aligned[kind][m['model_id']], aligned[kind][ref['id']]
                        columns = [i for i, (x, y) in enumerate(zip(a, b)) if x != '-' or y != '-']
                        filename = 'reference_alignments/'+gid+'/'+m['model_id']+'.'+ref['id']+'.'+kind+'.fasta'
                        target = root/filename; target.parent.mkdir(parents=True, exist_ok=True)
                        write_fasta(target, [(m['model_id']+' accession='+m['original_model']+' species='+m['species'].replace(' ', '_'), ''.join(a[i] for i in columns)),
                                             (ref['id']+' accession='+ref['accession']+' species='+ref['species'].replace(' ', '_'), ''.join(b[i] for i in columns))])
                        row['alignment_file'] = filename
                        if kind == 'protein' and float(pair['paired_positions']) > 0:
                            options.append((float(pair['identity'])*min(float(pair['coverage_model']), float(pair['coverage_reference'])), ref, row))
                    comparisons.append(row)
            best = max(options, key=lambda x: x[0]) if options else None
            m['closest_reference'] = dict(reference_id=best[1]['id'], species=best[1]['species'],
                accession=best[1]['accession'], symbol=best[1]['symbol'],
                identity=best[2]['identity'], coverage_model=best[2]['coverage_model'],
                coverage_reference=best[2]['coverage_reference']) if best else None
            sheets.append(dict(gene=group['name'], species=m['species'], copy_number=m['copy_number'] or '',
                model_id=m['model_id'], original_model=m['original_model'], source=m['source'],
                protein_aa=m['protein_aa'], sequence_guidance=m['guidance']['label'],
                guidance_detail=m['guidance']['detail'],
                closest_supplied_reference=best[1]['accession'] if best else '',
                reference_species=best[1]['species'] if best else '',
                protein_identity=best[2]['identity'] if best else '',
                protein_model_coverage=best[2]['coverage_model'] if best else '',
                protein_reference_coverage=best[2]['coverage_reference'] if best else '',
                protein_file=m.get('protein_file', ''), cds_file=m.get('cds_file', '')))
    fields = ['group','model_id','species','original_model','reference_id','reference_species','reference_accession','reference_symbol','kind','origin','self_reference','identity','coverage_model','coverage_reference','paired_positions','identical_positions','alignment_file']
    write_tsv(root/'reference_comparisons.tsv', comparisons, fields)
    write_tsv(root/'sequence_reference_sheet.tsv', sheets, list(sheets[0]) if sheets else ['gene','species','model_id'])
    from .stage3_handoff import all_alignments, model_decisions
    model_decisions(root,payload)
    from .stage3_finalization import prepare
    finalization=prepare(root,payload,models,references)
    summary = dict(payload, models=models, references=references, reference_comparisons=comparisons,
                   all_pairs=all_alignments(root,payload),final_choices=read_table(root/'final_model_choices.tsv'),finalization=finalization)
    serialized = json.dumps(summary, ensure_ascii=True, separators=(',', ':')).replace('<','\\u003c').replace('&','\\u0026')
    template = Path(__file__).with_name('stage3_summary.txt').read_text()
    final_script=Path(__file__).with_name('stage3_finalization.txt').read_text()
    (root/'report.html').write_text(template.replace('__DATA__',serialized).replace('__FINALIZATION_SCRIPT__',final_script))
