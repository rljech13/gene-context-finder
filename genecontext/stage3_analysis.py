"""Summaries of the exact reviewed sequences; no new identity assignments."""
from collections import Counter, defaultdict
import csv
from itertools import combinations
from pathlib import Path
from statistics import median

from .sequence import read_fasta, write_fasta, write_tsv
from .stage2_choices import constraints


def read_table(path):
    if not Path(path).exists(): return []
    with Path(path).open() as handle: return list(csv.DictReader(handle, delimiter='\t'))


def number(value):
    return float(value) if value not in ('', None) else None


def native_interval(row, start, end):
    return ((row['window_end']-end, row['window_end']-start) if row['strand']=='-'
            else (row['window_start']+start, row['window_start']+end))


def pair_relation(a, b):
    same=a['genome_id']==b['genome_id']
    overlap=max(0,min(a['window_end'],b['window_end'])-max(a['window_start'],b['window_start'])) if same and a['contig']==b['contig'] else 0
    return dict(pair_scope='within_genome' if same else 'between_genomes',
        genomic_window_relation='overlapping_native_windows' if overlap else 'distinct_native_windows',
        native_window_overlap_bp=overlap)


def extract_flanks(row, native_window, requested=50000):
    """Flanks start at the selected coding-model envelope, not a guessed TSS."""
    from .sequence import reverse_complement
    ws,we=row['window_start'],row['window_end'];start,end=row['start'],row['end']
    if start<ws or end>we: return {},dict(flank_status='model_not_fully_observed',upstream_bp=0,downstream_bp=0,flank_requested_bp=requested,upstream_native_start='',upstream_native_end='',downstream_native_start='',downstream_native_end='')
    low=native_window[max(ws,start-requested)-ws:start-ws]
    high=native_window[end-ws:min(we,end+requested)-ws]
    up,down=(reverse_complement(high),reverse_complement(low)) if row['strand']=='-' else (low,high)
    low_interval=(max(ws,start-requested),start);high_interval=(end,min(we,end+requested))
    upstream,downstream=(high_interval,low_interval) if row['strand']=='-' else (low_interval,high_interval)
    return dict(upstream=up,downstream=down),dict(flank_status='complete_requested_flanks' if len(up)==len(down)==requested else 'limited_by_saved_window',
        upstream_bp=len(up),downstream_bp=len(down),flank_requested_bp=requested,upstream_native_start=upstream[0],upstream_native_end=upstream[1],downstream_native_start=downstream[0],downstream_native_end=downstream[1])


def analyze(root, data, selected, catalogue, dna, blocks, categories, context):
    root=Path(root);by_id={r['locus_id']:r for r in selected};meta=constraints(data)
    physical=defaultdict(list)
    for m in data['models']:physical[meta[m['model_id']]['physical_id']].append(m)
    qc=[];pair_rows=[];paralogs=[];groups=[];profiles=[]
    for group in dict.fromkeys(r['group'] for r in selected):
        chosen=[r for r in selected if r['group']==group];folder=root/'queries'/group
        protein=read_table(folder/'protein_pairs.tsv');coding=read_table(folder/'cds_pairs.tsv')
        pmap={frozenset([r['sequence_a'],r['sequence_b']]):r for r in protein}
        cmap={frozenset([r['sequence_a'],r['sequence_b']]):r for r in coding}
        for a,b in combinations(chosen,2):
            p=pmap.get(frozenset([a['locus_id'],b['locus_id']]),{});c=cmap.get(frozenset([a['locus_id'],b['locus_id']]),{})
            row=dict(group=group,model_a=a['locus_id'],model_b=b['locus_id'],species_a=a['genome_id'],species_b=b['genome_id'],
                copy_a=a['copy_label'],copy_b=b['copy_label'],protein_identity=p.get('identity',''),protein_paired_positions=p.get('paired_positions',''),
                protein_coverage_a=p.get('coverage_a','') if p.get('sequence_a')==a['locus_id'] else p.get('coverage_b',''),
                protein_coverage_b=p.get('coverage_b','') if p.get('sequence_a')==a['locus_id'] else p.get('coverage_a',''),
                cds_identity=c.get('identity',''),cds_paired_positions=c.get('paired_positions',''),
                cds_coverage_a=c.get('coverage_a','') if c.get('sequence_a')==a['locus_id'] else c.get('coverage_b',''),
                cds_coverage_b=c.get('coverage_b','') if c.get('sequence_a')==a['locus_id'] else c.get('coverage_a',''),**pair_relation(a,b))
            pair_rows.append(row)
            if row['pair_scope']=='within_genome':paralogs.append(row)
        for m in chosen:
            alternatives=physical[m['physical_id']];opposite=[a for a in alternatives if a['source']!=m['source']]
            agree='exact_CDS_match' if any(a['cds_sha256']==m['cds_sha256'] for a in opposite) else 'different_CDS_models' if opposite else 'other_source_unavailable'
            ident=m['identification'];best=max(m['target_support'],key=lambda h:h['bitscore'],default={});flags=[]
            if ident['status']!='homology_assignment':flags.append(ident['status'])
            if agree=='different_CDS_models':flags.append('annotation_prediction_disagreement')
            if m['window_edge']:flags.append('window_edge')
            if m['coding_status']!='complete_predicted':flags.append(m['coding_status'])
            if 'coverage_query' in best and 'coverage_subject' in best and min(best['coverage_query'],best['coverage_subject'])<.5:flags.append('best_query_match_covers_less_than_half_of_one_sequence')
            candidates=[]
            for p in protein:
                if m['locus_id'] not in (p['sequence_a'],p['sequence_b']):continue
                other=p['sequence_b'] if p['sequence_a']==m['locus_id'] else p['sequence_a']
                if other.startswith('reference_') and number(p['identity']) is not None:
                    ca=number(p['coverage_a']) or 0;cb=number(p['coverage_b']) or 0
                    candidates.append((number(p['identity'])*min(ca,cb),other,p,ca if p['sequence_a']==m['locus_id'] else cb,cb if p['sequence_a']==m['locus_id'] else ca))
            ref=max(candidates,key=lambda x:x[0]) if candidates else None
            qc.append(dict(model_id=m['locus_id'],genome_id=m['genome_id'],copy_label=m['copy_label'],source=m['source'],original_model=m['original_model'],
                protein_aa=m['protein_length'],coding_status=m['coding_status'],cds_segments=len(m['cds_blocks']),
                annotation_name=m['annotation_name'],database_name=ident.get('proposed_name',''),database_status=ident['status'],source_agreement=agree,
                best_query=best.get('subject',''),query_identity=best.get('identity',''),query_model_coverage=best.get('coverage_query',''),query_reference_coverage=best.get('coverage_subject',''),
                best_aligned_reference=ref[1] if ref else '',aligned_reference_identity=ref[2]['identity'] if ref else '',aligned_model_coverage=ref[3] if ref else '',aligned_reference_coverage=ref[4] if ref else '',
                working_sequence_status='user_selected_with_cautions' if flags else 'user_selected_coding_consistent',cautions='; '.join(flags),native_DNA_CDS_protein_check='passed'))
        values=[number(p['protein_identity']) for p in pair_rows if p['group']==group and p['protein_identity']!='']
        high=[p for p in pair_rows if p['group']==group and min(number(p['protein_coverage_a']) or 0,number(p['protein_coverage_b']) or 0)>=.8]
        groups.append(dict(group=group,name=chosen[0]['group_name'],copies=len(chosen),species=len({m['genome_id'] for m in chosen}),expected_pairs=len(chosen)*(len(chosen)-1)//2,
            measured_pairs=len(values),protein_identity_min=min(values) if values else '',protein_identity_median=median(values) if values else '',protein_identity_max=max(values) if values else '',
            pairs_with_80_percent_coverage_both=len(high),protein_identity_median_high_coverage=median(number(p['protein_identity']) for p in high) if high else ''))
        alignment_path=folder/'proteins.aligned.faa'
        if alignment_path.exists():
            alignment=read_fasta(alignment_path);ids=[r['locus_id'] for r in chosen if r['locus_id'] in alignment];by_species=defaultdict(list)
            for m in chosen:
                if m['locus_id'] in alignment:by_species[m['genome_id']].append(m['locus_id'])
            for pos in range(len(next(iter(alignment.values())))):
                weights=Counter();observed=0
                for members in by_species.values():
                    aa=[alignment[mid][pos] for mid in members if alignment[mid][pos] in 'ACDEFGHIKLMNPQRSTVWY']
                    if aa:
                        observed+=1
                        for residue in aa:weights[residue]+=1/len(aa)
                residue,score=weights.most_common(1)[0] if weights else ('',0)
                profiles.append(dict(group=group,alignment_position=pos+1,observed_species=observed,eligible_species=len(by_species),species_coverage=observed/len(by_species) if by_species else '',consensus=residue,species_balanced_consensus_fraction=score/observed if observed>=2 else ''))
    for name,rows,fields in [('model_quality.tsv',qc,['model_id','cautions']),('target_pairs.tsv',pair_rows,['model_a','model_b']),('paralog_pairs.tsv',paralogs,['model_a','model_b']),('conservation_summary.tsv',groups,['group','copies']),('species_conservation_profile.tsv',profiles,['group','alignment_position'])]:
        write_tsv(root/name,rows,list(rows[0]) if rows else fields)
    for row in dna:
        a,b=by_id[row['locus_a']],by_id[row['locus_b']];row.update(pair_relation(a,b))
    for row in blocks:
        a,b=by_id[row['locus_a']],by_id[row['locus_b']]
        for side,m in [('a',a),('b',b)]:
            start,end=native_interval(m,row['start_'+side],row['end_'+side]);row['native_start_'+side]=start;row['native_end_'+side]=end
        same=a['genome_id']==b['genome_id'] and a['contig']==b['contig']
        overlap=max(0,min(row['native_end_a'],row['native_end_b'])-max(row['native_start_a'],row['native_start_b'])) if same else 0
        row['shared_native_sequence']=same and overlap/max(1,min(row['native_end_a']-row['native_start_a'],row['native_end_b']-row['native_start_b']))>=.9
    # Retain the original empty-table schemas if no comparisons were possible.
    if dna:write_tsv(root/'dna_pairs.tsv',dna,list(dna[0]))
    if blocks:write_tsv(root/'dna_blocks.tsv',blocks,list(blocks[0]))
    return dict(groups=groups,quality=qc,pairs=pair_rows,paralogs=paralogs,species_profiles=profiles,
        overlapping_window_pairs=sum(r['genomic_window_relation']=='overlapping_native_windows' for r in dna),
        reviewed_neighbor_models=len(context),reviewed_cross_locus_pairs=len(data.get('_reviewed_pairs',{})),
        source_agreement=dict(Counter(r['source_agreement'] for r in qc)),naming_status=dict(Counter(r['database_status'] for r in qc)))
