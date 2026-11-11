"""Compare explicitly reviewed locus-first models without discovery/prediction."""
from collections import defaultdict
from html import escape
import hashlib
import json
from pathlib import Path
import shutil
import time

from .alignment import backtranslate
from .comparison import orient_interval
from .comparison_report import protein_comparison, dna_svg, synteny_svg
from .context_compare import synteny_pairs
from .dna_compare import compare_dna
from .sequence import read_fasta, write_fasta, write_tsv, reverse_complement, translate, sha256
from .stage1 import dump, new_output, now
from .stage2_choices import verify_checkpoint, validate_choices, issue_text, constraints, copy_label
from .tools import Runner, tool_versions
from .stage3_analysis import analyze, extract_flanks
from .stage3_context import write_context


class TimedRunner(Runner):
    def run(self, args, output, log=None):
        started=time.perf_counter()
        try:return super().run(args,output,log)
        finally:
            if self.commands:self.commands[-1]['elapsed_seconds']=time.perf_counter()-started


def resource_usage():
    import resource,sys
    own=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
    divisor=1024 if sys.platform=='darwin' else 1
    return dict(process_peak_rss_kib=own.ru_maxrss/divisor,max_child_peak_rss_kib=children.ru_maxrss/divisor,
                note='Process high-water marks; not the aggregate cgroup peak. Host2 launcher records aggregate memory limits separately.')


def species(gid):
    return gid.split('_GC')[0].replace('_', ' ')


def matrix_html(path, selected):
    import csv
    if not path.exists():return '<p>Pairwise comparison unavailable; inspect errors.</p>'
    with path.open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
    pairs={frozenset([r['sequence_a'],r['sequence_b']]):r for r in rows}
    copies=defaultdict(int)
    for r in selected:copies[r['genome_id']]+=1
    out=['<div class="scroll"><table><tr><th>Species / copy</th>']
    for r in selected:
        words=species(r['genome_id']).split()
        label=(words[0][0]+'. '+words[1]) if len(words)>1 else words[0]
        copy=r.get('copy_label',r['locus_id'])
        full=species(r['genome_id'])+' · '+copy
        out.append('<th title="'+escape(full,quote=True)+'">'+escape(label)+(('<br><small>'+escape(copy)+'</small>') if copies[r['genome_id']]>1 else '')+'</th>')
    out.append('</tr>')
    for a in selected:
        out.append('<tr><th>'+escape(species(a['genome_id']))+(('<br><small>'+escape(a.get('copy_label',a['locus_id']))+'</small>') if copies[a['genome_id']]>1 else '')+'</th>')
        for b in selected:
            r=pairs.get(frozenset([a['locus_id'],b['locus_id']]))
            if a['locus_id']==b['locus_id']:out.append('<td>—</td>')
            elif r and r['identity']!='':
                tip=f"{r['paired_positions']} comparable positions; coverage {float(r['coverage_a'] or 0):.1%} / {float(r['coverage_b'] or 0):.1%} (table sequence order)"
                out.append('<td title="'+escape(tip,quote=True)+'">'+f"{float(r['identity']):.1%}"+'</td>')
            else:out.append('<td>Unavailable</td>')
        out.append('</tr>')
    return ''.join(out)+'</table></div>'


def prepare_inputs(data, state, source):
    """Validate exact saved sequence/model identity before any output is created."""
    problems = validate_choices(data, state)
    if problems:
        raise ValueError('STOP-2 decisions are not ready:\n' + '\n'.join(issue_text(data, p) for p in problems))
    proteins = read_fasta(source/'proteins.faa') if (source/'proteins.faa').stat().st_size else {}
    coding = read_fasta(source/'cds.fna') if (source/'cds.fna').stat().st_size else {}
    windows = read_fasta(source/'windows.fna', preserve_case=True)
    regions = {r['region_id']:r for r in data['regions']}
    metadata = constraints(data); selected = []
    for m in data['models']:
        c = state.get('models', {}).get(m['model_id'], {})
        if c.get('decision') != 'inspect' or c.get('role') != 'target':
            continue
        mid = m['model_id']; dna = coding.get(mid, ''); protein = proteins.get(mid, '')
        if not protein or '*' in protein or any(hashlib.sha256(seq.encode()).hexdigest() != m[field]
                for seq,field in [(dna,'cds_sha256'),(protein,'protein_sha256')]):
            raise ValueError('Selected CDS/protein hash or integrity mismatch: ' + mid)
        without_stop = dna[:-3] if translate(dna).endswith('*') else dna
        backtranslate({mid:protein},{mid:without_stop})
        r = regions[m['region_id']]; sequence = windows[r['region_id']]
        if len(sequence) != r['window_end']-r['window_start']:
            raise ValueError('Saved genomic window length differs: ' + r['region_id'])
        if all(r['window_start'] <= a < b <= r['window_end'] for a,b,_ in m['cds_blocks']):
            # The generic coding reader handles split codons and separate stops.
            from .sequence import Model, coding_sequence
            class SavedWindow:
                def fetch(self, contig, a, b, strand='+'):
                    text = sequence[a-r['window_start']:b-r['window_start']].upper()
                    if len(text) != b-a: raise ValueError('Selected coding segment outside saved window')
                    return reverse_complement(text) if strand == '-' else text
            model = Model(mid,m['contig'],m['strand'],cds=[tuple(x) for x in m['cds_blocks']],
                          stops=[tuple(x) for x in m['stop_blocks']])
            native, aa, _, _ = coding_sequence(model,SavedWindow())
            if (native,aa) != (dna,protein):
                raise ValueError('Selected model fails saved native-DNA roundtrip: ' + mid)
        group_name = c.get('comparison_group',r.get('gene','target'))
        group = 'group_' + hashlib.sha256(group_name.encode()).hexdigest()[:16]
        names = list(dict.fromkeys(h.get('gene_symbol','') for h in m['identification'].get('near_best',[]) if h.get('gene_symbol')))
        label = m['identification'].get('proposed_name') or m['annotation_name'] or ' / '.join(names) or 'Unassigned candidate'
        own_anchor=True  # Final comparisons use the chosen model envelope, not the discovery seed.
        row = dict(m,locus_id=mid,group=group,group_name=group_name,copy_label=copy_label(m,c),
            target_name=label,status=m['coding_status'],window_start=r['window_start'],window_end=r['window_end'],
            target_start=m['start'] if own_anchor else r['anchor_start'],
            target_end=m['end'] if own_anchor else r['anchor_end'],
            discovery_seed_start=r['anchor_start'],discovery_seed_end=r['anchor_end'],
            window_bp=len(sequence),window_fasta='windows/'+mid+'.fna',
            physical_id=metadata[mid]['physical_id'],recommendation_reason=c.get('notes') or state['regions'][r['region_id']]['notes'],
            target_override=c.get('target_override',False),selection_kind=c.get('selection_kind',''),source_region_id=r['region_id'])
        row.update(selection_status='explicit_model_choice',anchor_kind='reviewed_target',context_status='reviewed_neighbor_models')
        row['target_relative_start'],row['target_relative_end']=orient_interval(row,row['target_start'],row['target_end'])
        row['relative_cds']=[orient_interval(row,a,b) for a,b,_ in m['cds_blocks']]
        row['softmasked_bases']=sum(x.islower() for x in sequence)
        selected.append((row,protein,dna,sequence if m['strand']=='+' else reverse_complement(sequence)))
    return selected


def reviewed_neighbors(data,state,selected):
    """Only explicitly selected neighboring models become comparison anchors."""
    meta=constraints(data);models={m['model_id']:m for m in data['models']};parent={}
    def find(x):
        parent.setdefault(x,x)
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for pair in state.get('pairs',{}).values():
        a,b=[meta[pair[k]]['physical_id'] for k in ['model_a','model_b']]
        parent[find(b)]=find(a)
    components=defaultdict(set)
    for mid,info in meta.items():
        if info['physical_id'] in parent:components[find(info['physical_id'])].add((models[mid]['genome_id'],info['physical_id']))
    ambiguous={key for key,values in components.items() if len({gid for gid,p in values})<len(values)}
    context=[]
    for row in selected:
        seen=set()
        for mid,c in state.get('models',{}).items():
            m=models[mid]
            if m['region_id']!=row['source_region_id'] or c.get('decision')!='inspect' or c.get('role')!='neighbor':continue
            physical=meta[mid]['physical_id']
            if physical in seen:continue
            seen.add(physical);a,b=orient_interval(row,m['start'],m['end'])
            key=find(physical) if physical in parent else ''
            label=m['identification'].get('proposed_name') or m['annotation_name'] or 'Unassigned neighbor'
            context.append(dict(locus_id=row['locus_id'],group=row['group'],genome_id=m['genome_id'],model_id=mid,
                physical_id=physical,name=label,contig=m['contig'],start=m['start'],end=m['end'],strand=m['strand'],
                relative_start=a,relative_end=b,relative_strand='+' if m['strand']==row['strand'] else '-',
                evidence='provided_gff3' if m['source']=='annotation' else 'Tiberius_prediction',
                homology_group=key if key and key not in ambiguous else '',
                correspondence_status='ambiguous_component' if key in ambiguous else 'reviewed_pair' if key else 'not_selected_as_pair'))
    return context


def run(args):
    started=time.perf_counter();timings={};source,data=verify_checkpoint(args.results)
    state=json.loads(Path(args.decisions).read_text())
    inputs=prepare_inputs(data,state,source)
    # Stage-2 BLAST queries are extracted model proteins; its subjects are
    # the supplied target references. Neighbors must not enter the target MSA.
    raw_queries=source/'target_support'/'subjects.faa'
    references=read_fasta(raw_queries) if raw_queries.exists() and raw_queries.stat().st_size else {}
    if set(references)-set(data.get('query_members',{})):
        raise ValueError('Saved target references do not match the stage-2 query panel')
    if any(row['coding_status'] != 'complete_predicted' for row,_,_,_ in inputs) and not args.include_partial:
        raise ValueError('Selected partial models require explicit --include-partial; they are not silently omitted')
    if args.validate_only:
        print(json.dumps(dict(valid=True,review_id=data['review_id'],selected_models=len(inputs),
            selected_species=len({r['genome_id'] for r,_,_,_ in inputs}),output_created=False)))
        return 0
    versions=tool_versions(['mafft',args.dna_engine]) if len(inputs)>1 else {}
    if any(not v['path'] for v in versions.values()):
        raise ValueError('Missing comparison tool: '+', '.join(k for k,v in versions.items() if not v['path']))
    timings['input_validation_seconds']=time.perf_counter()-started
    order={r['region_id']:i for i,r in enumerate(data['regions'])}
    inputs.sort(key=lambda item:(order[item[0]['source_region_id']],item[0]['start'],item[0]['model_id']))
    root=new_output(args.out);runner=TimedRunner(root);errors=[];step=time.perf_counter()
    (root/'windows').mkdir();(root/'sequences').mkdir()
    shutil.copyfile(source/'evidence.json',root/'source.evidence.json')
    # Keep original evidence byte-for-byte and carry its separate support
    # supplements into the handoff; later ledgers must use the same constraints.
    if (source/'reconstructions.json').exists():
        shutil.copyfile(source/'reconstructions.json',root/'source.reconstructions.json')
    if (source/'annotation_quality.json').exists():
        shutil.copyfile(source/'annotation_quality.json',root/'source.annotation_quality.json')
    shutil.copyfile(args.decisions,root/'decisions.snapshot.json')
    selected=[];windows={};overrides={};catalogue=[];all_flanks=dict(upstream=[],downstream=[])
    for row,protein,cds,sequence in inputs:
        mid=row['model_id'];selected.append(row);windows[mid]=sequence
        overrides[mid]=dict(protein=protein,cds=cds)
        write_fasta(root/row['window_fasta'],[(mid,sequence)])
        for suffix,value in [('protein.faa',protein),('cds.fna',cds)]:
            write_fasta(root/'sequences'/(mid+'.'+suffix),[(mid,value)])
        native_window=sequence if row['strand']=='+' else reverse_complement(sequence)
        fully_observed=row['window_start']<=row['start']<row['end']<=row['window_end']
        gene=''
        if fully_observed:
            gene=native_window[row['start']-row['window_start']:row['end']-row['window_start']]
            if row['strand']=='-':gene=reverse_complement(gene)
            write_fasta(root/'sequences'/(mid+'.gene.fna'),[(mid,gene)])
        flanks,flank_info=extract_flanks(row,native_window)
        for kind,sequence_flank in flanks.items():
            write_fasta(root/'sequences'/(mid+'.'+kind+'.fna'),[(mid,sequence_flank)])
            all_flanks[kind].append((mid,sequence_flank))
        catalogue.append(dict(**flank_info,target=row['group_name'],species=species(row['genome_id']),genome_id=row['genome_id'],
            assembly_accession=('GC'+row['genome_id'].rsplit('_GC',1)[1]) if '_GC' in row['genome_id'] else '',
            physical_id=row['physical_id'],copy_label=row['copy_label'],selection_kind=row['selection_kind'],region_id=row['source_region_id'],model_id=mid,source=row['source'],
            original_model=row['original_model'],original_gene=row['original_gene'],gene_label=row['target_name'],
            contig=row['contig'],start=row['start'],end=row['end'],strand=row['strand'],protein_aa=len(protein),
            coding_status=row['coding_status'],window_edge=row['window_edge'],target_override=row['target_override'],
            protein_sha256=row['protein_sha256'],cds_sha256=row['cds_sha256'],
            genomic_window_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            gene_dna_status='exported' if fully_observed else 'not_fully_observed_in_saved_window',
            protein_file='sequences/'+mid+'.protein.faa',cds_file='sequences/'+mid+'.cds.fna',
            gene_file='sequences/'+mid+'.gene.fna' if fully_observed else '',window_file=row['window_fasta'],
            upstream_file='sequences/'+mid+'.upstream.fna' if 'upstream' in flanks else '',downstream_file='sequences/'+mid+'.downstream.fna' if 'downstream' in flanks else '',
            recommendation_reason=row['recommendation_reason'],comparison_scope='reviewed target candidates; copy orthology not inferred by this command'))
    for name,values in [('selected.proteins.faa',[(r['model_id'],p) for r,p,_,_ in inputs]),
                        ('selected.cds.fna',[(r['model_id'],c) for r,_,c,_ in inputs]),
                        ('selected.genomic_windows.fna',list(windows.items()))]:write_fasta(root/name,values)
    for kind,records in all_flanks.items():write_fasta(root/('selected.'+kind+'.fna'),records)
    fields=list(catalogue[0]) if catalogue else ['target','species','model_id','source','recommendation_reason']
    write_tsv(root/'sequence_catalogue.tsv',catalogue,fields)
    chosen_by_region={r['source_region_id'] for r in selected}
    observations=[]
    for r in data['regions']:
        decision=state['regions'][r['region_id']]
        observations.append(dict(genome_id=r['genome_id'],species=species(r['genome_id']),region_id=r['region_id'],
            decision=decision['decision'],compared=r['region_id'] in chosen_by_region,reason=decision['notes']))
    for gid in data.get('selection',{}).get('omitted_unsaved_genomes',[]):
        observations.append(dict(genome_id=gid,species=species(gid),region_id='',decision='stage1_review_pending',compared=False,reason='Not forwarded: stage 1 decision was not saved'))
    write_tsv(root/'review_summary.tsv',observations,['genome_id','species','region_id','decision','compared','reason'])
    queries=[dict(id=qid,protein=seq,coding='') for qid,seq in references.items()]
    groups=[dict(id=gid,name=next(r['group_name'] for r in selected if r['group']==gid),queries=list(references)) for gid in dict.fromkeys(r['group'] for r in selected)]
    timings['sequence_exports_seconds']=time.perf_counter()-step;step=time.perf_counter()
    protein_summary=protein_comparison(args,source,root,groups,queries,selected,runner,errors,sequence_overrides=overrides)
    timings['protein_and_codon_alignment_seconds']=time.perf_counter()-step;step=time.perf_counter()
    dna,blocks,categories=compare_dna(args,root,selected,windows,runner,errors)
    timings['genomic_alignment_seconds']=time.perf_counter()-step;step=time.perf_counter()
    context=reviewed_neighbors(data,state,selected);synteny=synteny_pairs(selected,context)
    write_tsv(root/'reviewed_neighbors.tsv',context,['locus_id','genome_id','model_id','physical_id','name','contig','start','end','strand','evidence','homology_group','correspondence_status'])
    write_tsv(root/'synteny.tsv',synteny,['group','locus_a','locus_b','unique_anchors','order_concordance','relative_orientation_concordance','status'])
    analysis=analyze(root,data,selected,catalogue,dna,blocks,categories,context)
    analysis['reviewed_cross_locus_pairs']=len(state.get('pairs',{}))
    write_context(root,source,data,state,selected)
    timings['context_and_summaries_seconds']=time.perf_counter()-step
    dump(root/'comparison.json',dict(groups=groups,selected=selected,protein=protein_summary,dna=dna,blocks=blocks,categories=categories,context=context,synteny=synteny,errors=errors,analysis=analysis))
    manifest=dict(stage='stage3_comparison',status='completed_with_errors' if errors else 'completed' if selected else 'no_selected_targets',
        created=now(),source_checkpoint=str(source),source_review_id=data['review_id'],
        source_evidence_sha256=sha256(source/'evidence.json'),decisions_sha256=sha256(root/'decisions.snapshot.json'),
        selected_models=len(selected),selected_species=len({r['genome_id'] for r in selected}),
        target_reference_ids=list(references),
        control_run=args.control_run,biological_approval_asserted=False,parameters=vars(args),tools=versions,
        commands=runner.commands,errors=errors,elapsed_seconds=time.perf_counter()-started,timings=timings,resources=resource_usage(),
        source_code_sha256={n:sha256(Path(__file__).with_name(n)) for n in ['stage3.py','stage3_analysis.py','stage3_report.py','stage3_report.txt','stage3_summary.py','stage3_summary.txt','stage3_handoff.py','stage3_finalization.py','stage3_finalization.txt','gene_identity.py','stage3_context.py','stage2_choices.py','alignment.py','comparison_report.py','dna_compare.py','stage2_neighborhood.txt','stage2_cluster.txt']})
    step=time.perf_counter()
    from .stage3_report import write_report as full_report
    full_report(root,data,catalogue,groups,selected,observations,blocks,context,errors,args.control_run,analysis,manifest)
    timings['report_seconds']=time.perf_counter()-step;manifest['elapsed_seconds']=time.perf_counter()-started;manifest['resources']=resource_usage()
    dump(root/'manifest.json',manifest)
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    print('Stage 3 report: '+str(root/'report.html'))
    return 1 if errors else 0


