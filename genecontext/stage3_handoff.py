"""Pairwise MSA projections and a faithful, read-only STOP-2 choice ledger."""
from collections import defaultdict
import json
from pathlib import Path

from .sequence import write_fasta, write_tsv
from .stage2_choices import constraints
from .stage3_analysis import read_table


def all_alignments(root, payload):
    root = Path(root)
    pairs = []
    labels={m['model_id']:(m['species'],m['original_model']) for m in payload['models']}
    for q,members in payload.get('query_members',{}).items():
        member=(members or [{}])[0]
        labels['reference_'+q]=(member.get('species',''),member.get('protein') or member.get('query') or q)
    for group in payload['groups']:
        gid = group['id']; folder = root/'queries'/gid
        for kind, table in [('protein', 'protein_pairs.tsv'), ('cds', 'cds_pairs.tsv')]:
            alignment = payload['alignments'].get(gid, {}).get(kind, {})
            for pair in read_table(folder/table):
                a, b = pair['sequence_a'], pair['sequence_b']
                if a not in alignment or b not in alignment:
                    raise ValueError('Pair table names a sequence absent from its alignment')
                sa, sb = alignment[a], alignment[b]
                columns = [i for i, (x, y) in enumerate(zip(sa, sb)) if x != '-' or y != '-']
                filename = f'all_alignments/{gid}/{kind}/{a}.{b}.fasta'
                file = root/filename; file.parent.mkdir(parents=True, exist_ok=True)
                headers={id:id+' accession='+labels.get(id,('',id))[1]+' species='+labels.get(id,('',id))[0].replace(' ','_') for id in (a,b)}
                write_fasta(file, [(headers[a], ''.join(sa[i] for i in columns)), (headers[b], ''.join(sb[i] for i in columns))])
                pairs.append(dict(group=gid, kind=kind, **pair, species_a=labels.get(a,('',''))[0],accession_a=labels.get(a,('',''))[1],
                                  species_b=labels.get(b,('',''))[0],accession_b=labels.get(b,('',''))[1],alignment_file=filename))
    fields = ['group','kind','sequence_a','sequence_b','species_a','accession_a','species_b','accession_b','paired_positions','identical_positions','identity','coverage_a','coverage_b','alignment_file']
    write_tsv(root/'all_alignment_pairs.tsv', pairs, fields)
    return pairs


def model_decisions(root, payload):
    """Keep outcomes literal; an unchosen alternative is not a biological rejection."""
    root = Path(root)
    data = json.loads((root/'source.evidence.json').read_text())
    if (root/'source.annotation_quality.json').exists():
        from .stage2_quality import with_annotation_quality
        data = with_annotation_quality(root,data,evidence_name='source.evidence.json',
                                       supplement_name='source.annotation_quality.json')
    if (root/'source.reconstructions.json').exists():
        from .stage2_reconstruct import load as load_reconstructions
        # The loader binds metadata to evidence.json; use a read-only locator.
        data = load_reconstructions(root, data, evidence_name='source.evidence.json',
                                   supplement_name='source.reconstructions.json')
    state = json.loads((root/'decisions.snapshot.json').read_text())
    metadata = constraints(data)
    regions = {r['region_id']: r for r in data['regions']}
    chosen = {m['model_id']: m for m in payload['models']}
    group_by_name = {g['name']: g['id'] for g in payload['groups']}
    physical_chosen = defaultdict(list)
    seeds = set()
    for m in data['models']:
        mid = m['model_id']; choice = state.get('models', {}).get(mid, {})
        physical = metadata[mid]['physical_id']
        if (m.get('overlaps_anchor') and m.get('same_strand_as_anchor', True)) or mid in chosen or choice.get('decision') in {'reject','not_selected','unresolved'}:
            seeds.add(physical)
        if mid in chosen:
            physical_chosen[physical].append(chosen[mid])
    observations = []
    for m in data['models']:
        mid=m['model_id']; rid=m['region_id']; region=regions[rid]
        choice=state.get('models', {}).get(mid, {}); rc=state.get('regions', {}).get(rid, {})
        physical=metadata[mid]['physical_id']; targets=physical_chosen[physical]
        group_name=choice.get('comparison_group') or region.get('gene', 'target')
        if targets: group_name=next(g['name'] for g in payload['groups'] if g['id']==targets[0]['group'])
        group=group_by_name.get(group_name, '')
        note=choice.get('notes','').strip(); region_note=rc.get('notes','').strip()
        if mid in chosen:
            status='Retained'; reason=note or region_note or 'Selected as target in the saved STOP-2 review.'
        elif choice.get('decision')=='inspect' and choice.get('role')=='neighbor':
            status='Context only'; reason=note or 'Selected as a neighboring gene, not as a target sequence.'
        elif physical not in seeds:
            status='Context only'; reason='Available neighboring model; not a target candidate in this review.'
        elif targets:
            status='Not retained'
            identical=[c for c in targets if c['cds_sha256']==m['cds_sha256']]
            reason=note or (('Identical CDS retained as ' if identical else 'Alternative model; selected model: ')+', '.join(c['original_model'] for c in (identical or targets))+'.')
        elif rc.get('decision')=='unresolved' or choice.get('decision')=='unresolved':
            status='Unresolved'; reason=note or region_note or 'No target model was resolved at STOP 2.'
        elif rc.get('decision')=='exclude':
            status='Not retained'; reason='Locus excluded at STOP 2.'+(' '+(note or region_note) if note or region_note else '')
        else:
            status='Not retained'; reason=note or 'Not selected as a target in the saved STOP-2 decisions.'
        species=m['genome_id'].split('_GC')[0].replace('_',' ')
        observations.append(dict(group=group, gene=group_name, species=species, genome_id=m['genome_id'],
            physical_id=physical, region_id=rid, model_id=mid, original_model=m['original_model'],
            source=m['source'], gene_label=m['identification'].get('proposed_name') or m.get('annotation_name',''),
            contig=m['contig'], start=m['start'], end=m['end'], strand=m['strand'],
            protein_aa=m['protein_length'], cds_sha256=m['cds_sha256'], status=status, reason=reason,
            saved_model_decision=choice.get('decision',''), saved_model_role=choice.get('role',''),
            saved_model_notes=note, saved_region_decision=rc.get('decision',''), saved_region_notes=region_note,
            automatic_model_comments='; '.join(v for v in state.get('automatic_comments',{}).get(rid,[]) if isinstance(v,str)),
            cds_file=chosen[mid].get('cds_file','') if mid in chosen else '',
            protein_file=chosen[mid].get('protein_file','') if mid in chosen else ''))
    fields=list(observations[0]) if observations else ['model_id','status','reason']
    write_tsv(root/'model_decision_audit.tsv', observations, fields)
    # One compact row per physical target candidate. Identical observations in
    # overlapping windows remain in the audit, not as duplicate final choices.
    buckets=defaultdict(list)
    for row in observations:
        if row['physical_id'] in seeds and row['status']!='Context only':
            buckets[row['group'],row['physical_id']].append(row)
    final=[]; summary=[]
    for (group, physical), rows in buckets.items():
        by_model=defaultdict(list)
        for row in rows:by_model[row['source'],row['original_model'],row['cds_sha256']].append(row)
        alternatives=[]
        for repeated in by_model.values():
            row=next((r for r in repeated if r['status']=='Retained'),
                     next((r for r in repeated if r['status']=='Unresolved'), repeated[0]))
            alternatives.append(dict(row,observations=len(repeated)))
        retained=[r for r in alternatives if r['status']=='Retained']
        others=[r for r in alternatives if r['status']!='Retained']
        status='Retained' if retained else 'Unresolved' if any(r['status']=='Unresolved' for r in others) else 'Not retained'
        representative=(retained or [r for r in alternatives if r['status']=='Unresolved'] or alternatives)[0]
        row=dict(group=group,physical_id=physical,species=representative['species'],
                 contig=representative['contig'],start=min(r['start'] for r in rows),end=max(r['end'] for r in rows),
                 strand=representative['strand'],status=status,retained=retained,alternatives=others,
                 reason=representative['reason'],automatic_model_comments=representative['automatic_model_comments'],observations=len(rows))
        final.append(row)
        summary.append(dict(group=group,species=row['species'],contig=row['contig'],start=row['start'],end=row['end'],strand=row['strand'],status=status,
            retained_models='; '.join(r['original_model']+' ['+r['source']+']' for r in retained),
            other_models='; '.join(r['original_model']+' ['+r['source']+'] ('+r['status']+')' for r in others),reason=row['reason']))
    species_order={}
    for region in data['regions']:
        species=region['genome_id'].split('_GC')[0].replace('_',' ')
        species_order.setdefault(species,len(species_order))
    ordering=lambda r:(species_order.get(r['species'],len(species_order)),r['contig'],r['start'])
    final.sort(key=ordering);summary.sort(key=ordering)
    write_tsv(root/'model_choice_history.tsv', summary, list(summary[0]) if summary else ['species','status','retained_models','other_models'])
    # The handoff is the selected set only. Past alternatives and unresolved
    # loci belong to the separate history, not to the final working catalogue.
    retained_rows=[]
    automatic_by_model={r['model_id']:r['automatic_model_comments'] for r in observations}
    for model in payload['models']:
        retained_rows.append(dict(group=model['group'],gene=next(g['name'] for g in payload['groups'] if g['id']==model['group']),species=model['species'],
            model_id=model['model_id'],original_model=model['original_model'],
            copy_label=model.get('copy_label',''),source=model['source'],
            automatic_model_comments=automatic_by_model.get(model['model_id'],''),
            contig=model['contig'],start=model['start'],end=model['end'],strand=model['strand'],
            protein_aa=model.get('protein_aa',model.get('protein_length','')),
            cds_file=model.get('cds_file',''),protein_file=model.get('protein_file','')))
    write_tsv(root/'final_model_choices.tsv',retained_rows,list(retained_rows[0]) if retained_rows else ['model_id','species','original_model','cds_file','protein_file'])
    return final
