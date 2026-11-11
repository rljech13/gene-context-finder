"""Selected-locus reannotation and evidence, ending at human review STOP 2."""
from collections import defaultdict
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time
from uuid import uuid4

from . import __version__
from .comparison import table
from .context_predict import prediction_model, validate_prediction_inputs
from .early_synteny import datasets, parse_annotation, summarize_hsps
from .sequence import coding_sequence, open_text, reverse_complement, sha256, write_fasta, write_tsv
from .stage1 import dump, new_output, now, validate_decisions, verify
from .stage2_db import verify_db
from .tools import Runner, tool_versions


class TimedRunner(Runner):
    def run(self,args,output,log=None):
        started=time.perf_counter();before=len(self.commands)
        try:return super().run(args,output,log)
        finally:
            if len(self.commands)>before:self.commands[-1]['elapsed_seconds']=time.perf_counter()-started
            dump(self.root/'progress.json',dict(updated=now(),commands_completed=len(self.commands),
                last_command=self.commands[-1] if self.commands else None))

def ident(*values):
    return hashlib.sha256(json.dumps(values,sort_keys=True).encode()).hexdigest()[:24]


def selected_regions(data, args):
    """Validate existing decisions; never turn a draft checkbox into approval."""
    provenance = {}
    if args.selection:
        root = Path(args.selection).resolve()
        manifest = json.loads((root/'manifest.json').read_text())
        if manifest.get('review_id') != data['review_id'] or manifest.get('stage')!='approved_loci':
            raise ValueError('Selection belongs to another stage 1 checkpoint')
        path = root/'decisions.snapshot.tsv'
        if sha256(path)!=manifest['decisions_sha256']:
            raise ValueError('Selected decision snapshot changed')
        args.decisions = str(path)
        provenance['selection'] = str(root)
    omitted=[]; subset=data
    if args.draft:
        path=Path(args.draft).resolve(); draft=json.loads(path.read_text())
        if draft.get('schema')!='stage1_pair_draft_v1' or draft.get('review_id')!=data['review_id']:
            raise ValueError('Wrong stage 1 draft')
        pairs=draft.get('pairs',[]); genomes={g['genome_id'] for g in data['genomes']}
        if not isinstance(pairs,list) or any(not isinstance(p,dict) or type(p.get('committed')) is not bool for p in pairs):
            raise ValueError('Draft pairs require explicit boolean committed states')
        if len(pairs)!=len(genomes) or {p['genome_id'] for p in pairs}!=genomes:
            raise ValueError('Draft must contain each stage 1 genome exactly once')
        omitted=[p['genome_id'] for p in pairs if not p.get('committed')]
        if omitted and not args.committed_only:
            raise ValueError('Draft has unsaved pairs. Finish the review, or explicitly use --committed-only.')
        accepted={p['genome_id']:p for p in pairs if p.get('committed')}
        subset=dict(data,genomes=[g for g in data['genomes'] if g['genome_id'] in accepted],
                    loci=[l for l in data['loci'] if l['genome_id'] in accepted])
        rows=[]
        loci={l['locus_id']:l for l in data['loci']}
        for gid,p in accepted.items():
            selected=p.get('selected',[])
            if len(set(selected))!=len(selected) or any(loci.get(i,{}).get('genome_id')!=gid for i in selected):
                raise ValueError('Unknown or duplicate selection in draft: '+gid)
            if p.get('action')!='inspect' and selected:
                raise ValueError('Excluded/unresolved pair contains selections: '+gid)
            own=[l for l in subset['loci'] if l['genome_id']==gid]
            for l in own or [None]:
                chosen=l and l['locus_id'] in selected
                rows.append(dict(review_id=data['review_id'],genome_id=gid,locus_id=l['locus_id'] if l else '',
                    decision=('inspect' if chosen else 'exclude') if p['action']=='inspect' else p['action'],
                    pair_decision=p['action'],reason=p.get('reason',''),notes=p.get('notes',''),
                    **(p.get('overrides',{}).get(l['locus_id'],{}) if l else {})))
        provenance['draft'] = str(path)
    else:
        path=Path(args.decisions).resolve(); rows=table(path)
    result=validate_decisions(subset,rows)
    if not result['valid']:
        raise ValueError('Stage 1 decisions are not valid:\n'+'\n'.join(result['errors']))
    if not result['regions']:
        raise ValueError('No explicitly selected region')
    return result['regions'],dict(provenance,decision_source_sha256=sha256(path),omitted_unsaved_genomes=omitted,
        neighbor_shortlist=(draft.get('neighbor_shortlist',[]) if args.draft else [])),rows


def planned_regions(source, data, args):
    """One boundary policy for both dry planning and actual DNA extraction."""
    if args.committed_only and not args.draft:raise ValueError('--committed-only requires --draft')
    regions,selection,rows=selected_regions(data,args)
    if args.selection:
        original=json.loads((Path(args.selection)/'manifest.json').read_text())
        if original['source_evidence_sha256']!=sha256(source/'evidence.json'):
            raise ValueError('Selection source evidence differs')
    if args.only_genomes:
        if set(args.only_genomes)-{r['genome_id'] for r in regions}:
            raise ValueError('--only-genomes contains an unapproved/unknown genome')
        regions=[r for r in regions if r['genome_id'] in set(args.only_genomes)]
        selection['explicit_genome_subset']=args.only_genomes
    from .stage2_windows import resolve_windows
    from .stage2_intervals import canonical_windows
    policy=getattr(args,'window_policy','pre-synteny')
    regions=resolve_windows(data,regions,rows,policy)
    if any(r['window_end']-r['window_start']>args.max_window_bases for r in regions):
        raise ValueError('Resolved window exceeds --max-window-bases; review its bounds or explicitly raise the cap')
    intervals=canonical_windows(regions,args.max_window_bases)
    selection['window_policy']=policy
    return regions,selection,rows,intervals


def window_plan(source,data,args,regions,selection,intervals):
    return dict(schema='stage2_window_plan_v1',source_checkpoint=str(source),
        source_review_id=data['review_id'],source_evidence_sha256=sha256(source/'evidence.json'),
        selection=selection,max_window_bases=args.max_window_bases,regions=regions,
        prediction_intervals=intervals,total_region_bases=sum(r['window_end']-r['window_start'] for r in regions),
        total_prediction_bases=sum(r['window_end']-r['window_start'] for r in intervals),
        note='Planned bounds; extraction, prediction and human model approval are separate operations.')


def plan(args):
    source,data=verify(args.results)
    regions,selection,_,intervals=planned_regions(source,data,args)
    payload=window_plan(source,data,args,regions,selection,intervals)
    root=new_output(args.out)
    dump(root/'window_plan.json',payload)
    write_tsv(root/'windows.tsv',[dict(r,window_bp=r['window_end']-r['window_start'],
        preliminary_neighbors=len(r['window_resolution']['neighbor_genes']),
        window_reason=r['window_resolution']['reason']) for r in regions],
        ['genome_id','locus_id','region_id','contig','anchor_start','anchor_end','window_start','window_end',
         'window_bp','window_source','preliminary_neighbors','window_reason'])
    print('Window plan — '+str(root/'windows.tsv'))
    return 0


def extract_intervals(path, requests):
    """One streaming pass over ordinary/gzipped FASTA; no whole-genome expansion."""
    by_contig=defaultdict(list)
    for key,(contig,start,end) in requests.items():
        if not 0<=start<end:raise ValueError('Invalid extraction interval: '+key)
        by_contig[contig].append((key,start,end))
    pieces={key:[] for key in requests}; lengths={}; seen=set(); contig=None; offset=0; active=[]; cursor=0; pending=[]
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if contig is not None:lengths[contig]=offset
                words=line[1:].split()
                if not words:raise ValueError('Empty FASTA header')
                contig=words[0]
                if contig in seen:raise ValueError('Duplicate FASTA contig: '+contig)
                seen.add(contig);offset=0;active=[];cursor=0;pending=sorted(by_contig.get(contig,[]),key=lambda x:x[1])
            else:
                seq=''.join(line.split())
                if seq and contig is None:raise ValueError('FASTA sequence before header')
                while cursor<len(pending) and pending[cursor][1]<offset+len(seq):
                    active.append(pending[cursor]);cursor+=1
                active=[item for item in active if item[2]>offset]
                for key,a,b in active:
                    lo=max(a,offset);hi=min(b,offset+len(seq))
                    if lo<hi:pieces[key].append(seq[lo-offset:hi-offset])
                offset+=len(seq)
    if contig is not None:lengths[contig]=offset
    result={key:''.join(parts) for key,parts in pieces.items()}
    for key,(contig,a,b) in requests.items():
        if contig not in lengths or b>lengths[contig] or len(result[key])!=b-a:
            raise ValueError('Missing contig or interval outside genome: '+key)
        if not set(result[key].upper())<=set('ACGTRYMKSWBDHVN'):
            raise ValueError('Non-DNA symbols in selected region: '+key)
    return result,lengths


class Sequences:
    def __init__(self, parts, requests):
        self.parts=parts;self.requests=requests
        self.exact={(c,a,b):parts[k] for k,(c,a,b) in requests.items()}
    def fetch(self,c,a,b,strand='+'):
        seq=self.exact[c,a,b].upper()
        return reverse_complement(seq) if strand=='-' else seq


class WindowSequence:
    def __init__(self,windows):self.windows=windows
    def fetch(self,c,a,b,strand='+'):
        seq=self.windows[c][a:b].upper()
        if len(seq)!=b-a:raise ValueError('CDS outside prediction window')
        return reverse_complement(seq) if strand=='-' else seq


def model_record(model, region, source, cds, protein, status, issues, gene_name=''):
    offset=region['window_start'] if source!='annotation' else 0
    start=min([model.start]+[a for a,b in model.stops])+offset
    end=max([model.end]+[b for a,b in model.stops])+offset
    group=ident(region['genome_id'],region['contig'],model.strand,
        ('annotation',model.gene) if source=='annotation' else ('prediction',[(a+offset,b+offset) for a,b,_ in model.cds]))
    return dict(model_id='M'+ident(region['region_id'],source,model.id),region_id=region['region_id'],
        locus_id=region['locus_id'],genome_id=region['genome_id'],gene_group=group,
        original_model=model.id,original_gene=model.gene,source=source,annotation_name=gene_name,
        contig=region['contig'],strand=model.strand,start=start,end=end,
        cds_blocks=[[a+offset,b+offset,p] for a,b,p in model.cds],
        stop_blocks=[[a+offset,b+offset] for a,b in model.stops],
        coding_status=status,issues=issues,protein_length=len(protein),
        annotation_flags=list(model.attrs.get('_quality_flags',[])) if source=='annotation' else [],
        cds_sha256=hashlib.sha256(cds.encode()).hexdigest(),protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
        overlaps_anchor=max(start,region['anchor_start'])<min(end,region['anchor_end']),
        same_strand_as_anchor=model.strand==region['strand'],
        window_edge=start<=region['window_start'] or end>=region['window_end'])


def blast_pairs(runner, root, queries, subjects, threads, evalue):
    if not queries or not subjects:return []
    root.mkdir(parents=True,exist_ok=True)
    write_fasta(root/'queries.faa',sorted(queries.items()));write_fasta(root/'subjects.faa',sorted(subjects.items()))
    runner.run(['makeblastdb','-in',root/'subjects.faa','-dbtype','prot','-out',root/'db'],root/'makeblastdb.stdout')
    runner.run(['blastp','-query',root/'queries.faa','-db',root/'db','-num_threads',threads,
        '-evalue',evalue,'-max_target_seqs',len(subjects),'-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'],root/'blast.tsv')
    return parse_hits(root/'blast.tsv',queries,{k:len(v) for k,v in subjects.items()})


def parse_hits(path, queries, subjects):
    grouped=defaultdict(list)
    with open(path) as handle:
        for line in handle:
            f=line.split()
            if len(f)!=12 or f[0] not in queries or f[1] not in subjects:raise ValueError('Unknown BLAST protein ID')
            qa,qb,sa,sb=int(f[4])-1,int(f[5]),int(f[6])-1,int(f[7])
            ql=len(queries[f[0]]);sl=subjects[f[1]]
            if not (0<=qa<qb<=ql and 0<=sa<sb<=sl and (int(f[10]),int(f[11]))==(ql,sl)):
                raise ValueError('BLAST protein coordinates/lengths differ')
            identity,evalue,bitscore=float(f[2])/100,float(f[8]),float(f[9])
            if (not all(math.isfinite(x) for x in [identity,evalue,bitscore]) or
                    not 0<=identity<=1 or evalue<0 or bitscore<0 or int(f[3])<=0):
                raise ValueError('Invalid BLAST protein statistics')
            grouped[f[0],f[1]].append(dict(qstart=qa,qend=qb,sstart=sa,send=sb,qlen=ql,slen=sl,
                length=int(f[3]),identity=identity,evalue=evalue,bitscore=bitscore))
    return [dict(query=q,subject=s,**summarize_hsps(h)) for (q,s),h in sorted(grouped.items())]


def identify(runner, root, proteins, db, info, args):
    root.mkdir()
    if not proteins:return {},[]
    write_fasta(root/'queries.faa',sorted(proteins.items()))
    runner.run(['blastp','-query',root/'queries.faa','-db',db/'blastdb','-num_threads',args.threads,
        '-evalue',args.identification_evalue,'-max_target_seqs',info['proteins'],
        '-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'],root/'blast.tsv')
    con=sqlite3.connect('file:'+str(db/'metadata.sqlite')+'?mode=ro',uri=True)
    lengths=dict(con.execute('SELECT id,length FROM proteins'))
    hits=parse_hits(root/'blast.tsv',proteins,lengths);accepted=defaultdict(list)
    for h in hits:
        if min(h['coverage_query'],h['coverage_subject'])>=args.identification_coverage and h['identity']>=args.identification_identity:
            accepted[h['query']].append(h)
    result={}
    for pid in proteins:
        candidates=sorted(accepted[pid],key=lambda h:(-h['bitscore'],h['evalue'],h['subject']))
        near=[h for h in candidates if h['bitscore']>=candidates[0]['bitscore']*args.competing_ratio] if candidates else []
        names=set();options=[]
        for h in near:
            rid,original,gene,desc,taxon,length=con.execute('SELECT * FROM proteins WHERE id=?',(h['subject'],)).fetchone()
            # Case conventions vary by taxon (AICDA/Aicda); paralog suffixes
            # and accession identifiers remain distinct.
            names.add(('symbol',gene.casefold()) if gene else ('accession',original))
            options.append(dict(h,accession=original,gene_symbol=gene,description=desc,taxon_id=taxon))
        result[pid]=dict(status='no_passing_match' if not near else 'ambiguous_family_assignment' if len(names)>1 else 'homology_assignment',
            proposed_name=(options[0]['gene_symbol'] or options[0]['accession']) if len(names)==1 else '',
            symbol_variants=sorted({h['gene_symbol'] for h in options if h['gene_symbol']}),
            name_equivalence='symbol_case_only; accessions_exact',passing_hits=len(candidates),near_best=options)
    con.close()
    return result,hits


def reidentify(args):
    """Replay only protein identification from a sealed STOP-2 checkpoint."""
    import shutil
    from .stage2_choices import verify_checkpoint
    source,data=verify_checkpoint(args.results);seal=json.loads((source/'seal.json').read_text())
    context_file=source/'annotation_context.json'
    if context_file.exists():
        context=json.loads(context_file.read_text())
        if context.get('schema')!='stage2_annotation_context_v1' or context.get('source_evidence_sha256')!=sha256(source/'evidence.json'):
            raise ValueError('Annotation context belongs to another scientific checkpoint')
        data['annotation_features']=context['features']
        data['provenance']['annotation_context']=dict(source_sha256=sha256(context_file),sources=context['sources'])
    db,info=verify_db(args.protein_db)
    target=Path(args.out).resolve()
    if target==source or source in target.parents or target in source.parents:raise ValueError('Use a separate result directory')
    for key in ['identification_coverage','identification_identity','identification_evalue','competing_ratio']:
        setattr(args,key,data['provenance']['parameters'][key])
    root=new_output(target)
    for name in seal:
        if Path(name).parts[0]=='identification' or name in {'identification_hits.json','report.html','presentation.json','annotation_context.json','annotation_quality.json','reconstructions.json','prior_synteny_context.json'}:continue
        path=root/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,path)
    from .sequence import read_fasta
    proteins=read_fasta(root/'proteins.faa');runner=TimedRunner(root);started=time.perf_counter()
    identified,hits=identify(runner,root/'identification',proteins,db,info,args)
    for m in data['models']:
        m['identification']=identified.get(m['model_id'],dict(status='no_usable_protein',proposed_name='',near_best=[]))
    data['provenance']['reidentification']=dict(source=str(source),source_evidence_sha256=sha256(source/'evidence.json'),
        source_review_id=data['review_id'],elapsed_seconds=time.perf_counter()-started,commands=runner.commands,
        implementation_sha256=sha256(Path(__file__)),note='Only broad-database identification repeated; genomic extraction, prediction, target-query and neighbor searches reused unchanged.')
    data['review_id']=str(uuid4());data['created']=now()
    data['errors']=[e for e in data['errors'] if e.get('stage')!='identification']
    data['status']='completed_with_errors' if data['errors'] else 'awaiting_manual_review'
    data['provenance']['protein_db']=info;data['provenance']['protein_db_path']=str(db)
    dump(root/'evidence.json',data)
    # Identification changes no coordinates or sequence hashes. Rebind validated
    # support records to the new evidence snapshot rather than copying stale IDs.
    for name in ['annotation_quality.json','reconstructions.json','prior_synteny_context.json']:
        if (source/name).exists():
            supplement=json.loads((source/name).read_text())
            key='stage2_evidence_sha256' if name=='prior_synteny_context.json' else 'source_evidence_sha256'
            supplement[key]=sha256(root/'evidence.json');dump(root/name,supplement)
    dump(root/'identification_hits.json',hits);write_outputs(root,data)
    print('STOP 2 — '+str(root/'report.html'));return 1 if data['errors'] else 0


def assign_groups(models):
    # A predicted transcript may be an alternative of ONE annotated gene.
    # Fusions overlapping several annotated genes stay unresolved, never bridge
    # those distinct genes into a single synteny anchor.
    from .sequence import union_length
    annotation=defaultdict(list)
    for m in models:
        if m['source']=='annotation':annotation[m['genome_id'],m['contig'],m['strand']].append(m)
    for m in models:
        if m['source']=='annotation':continue
        overlap=set()
        for a in annotation[m['genome_id'],m['contig'],m['strand']]:
            blocks=[(max(x[0],y[0]),min(x[1],y[1])) for x in m['cds_blocks'] for y in a['cds_blocks'] if max(x[0],y[0])<min(x[1],y[1])]
            if union_length(blocks)>=.5*min(union_length([b[:2] for b in m['cds_blocks']]),union_length([b[:2] for b in a['cds_blocks']])):
                overlap.add(a['gene_group'])
        m['overlapping_annotation_groups']=sorted(overlap)
        if len(overlap)==1:m['gene_group']=next(iter(overlap))
        elif len(overlap)>1:m['issues']=sorted(set(m['issues'])|{'overlaps_multiple_annotated_genes'})

def neighbor_pairs(models, hits, args):
    lookup={m['model_id']:m for m in models};best={}
    target_groups={m['gene_group'] for m in models if m.get('overlaps_anchor')}
    for h in hits:
        a,b=lookup[h['query']],lookup[h['subject']]
        if a['genome_id']==b['genome_id'] or min(h['coverage_query'],h['coverage_subject'])<args.identification_coverage or h['identity']<args.identification_identity:continue
        key=(a['gene_group'],b['gene_group'])
        if key not in best or h['bitscore']>best[key]['bitscore']:
            best[key]=dict(h,gene_a=a['gene_group'],gene_b=b['gene_group'],model_a=a['model_id'],model_b=b['model_id'])
    competing=defaultdict(list)
    for h in best.values():
        competing[h['gene_a'],lookup[h['model_b']]['genome_id']].append(h)
    out=[]
    for h in best.values():
        peers=competing[h['gene_a'],lookup[h['model_b']]['genome_id']]
        near=[p for p in peers if p['bitscore']>=max(x['bitscore'] for x in peers)*args.competing_ratio]
        if h not in near:continue
        reverse=competing.get((h['gene_b'],lookup[h['model_a']]['genome_id']),[])
        reciprocal=bool(reverse) and h['gene_a'] in {p['gene_b'] for p in reverse if p['bitscore']>=max(x['bitscore'] for x in reverse)*args.competing_ratio}
        out.append(dict(h,kind='target_similarity' if h['gene_a'] in target_groups or h['gene_b'] in target_groups else 'neighbor_proposal',status='unique_reciprocal_proposal' if len(near)==1 and reciprocal and sum(p['bitscore']>=max(x['bitscore'] for x in reverse)*args.competing_ratio for p in reverse)==1 else 'ambiguous_proposal',orthology_confirmed=False))
    return out



def extract_genome(gid, own, entry, prediction_intervals, lmap, source_genomes, digest):
    receipt=dict(genome_id=gid)
    fasta=Path(entry['fasta']);before=fasta.stat()
    known=source_genomes[gid].get('fasta_sha256')
    if known and known!=digest:raise ValueError('FASTA differs from stage 1 discovery input')
    requests={r['region_id']:(r['contig'],r['window_start'],r['window_end']) for r in own}
    for c in prediction_intervals:
        if c['genome_id']==gid:
            requests[c['region_id']]=(c['contig'],c['window_start'],c['window_end'])
    gene_names={};annotation=[];staged=[];staged_features=[]
    if entry.get('gff3'):
        genes,parsed,_=parse_annotation(entry['gff3']);gene_names={g:rec['name'] for g,rec in genes.items()}
        for r in own:
            for item in parsed.values():
                m=item['model']
                if m.contig==r['contig'] and max(m.start,r['window_start'])<min(m.end,r['window_end']):
                    annotation.append((r,m))
                    for a,b,_ in m.cds:requests['C'+ident(m.contig,a,b)]=(m.contig,a,b)
                    for a,b in m.stops:requests['C'+ident(m.contig,a,b)]=(m.contig,a,b)
        receipt['annotation_sha256']=sha256(entry['gff3'])
        expected=source_genomes[gid].get('annotation_inputs',{}).get('gff3',{}).get('sha256')
        if expected and expected!=receipt['annotation_sha256']:raise ValueError('Annotation differs from stage 1 dataset; use its matching release')
        from .annotation_context import read_features
        staged_features=read_features(entry['gff3'],own)
        del genes,parsed
    pieces,lengths=extract_intervals(fasta,requests)
    if (fasta.stat().st_size,fasta.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):raise ValueError('FASTA changed during extraction')
    extracted={}
    for r in own:
        if lengths[r['contig']]!=lmap[r['locus_id']]['contig_length']:raise ValueError('Contig length differs from stage 1')
        extracted[r['region_id']]=pieces[r['region_id']]
    for c in prediction_intervals:
        if c['genome_id']==gid:extracted[c['region_id']]=pieces[c['region_id']]
    genome=Sequences(pieces,requests)
    for r,m in annotation:
        cds,protein,status,issues=coding_sequence(m,genome)
        record=model_record(m,r,'annotation',cds,protein,status,issues,gene_names.get(m.gene,''))
        staged.append((record,cds,protein))
    receipt.update(status='completed',fasta_sha256=digest,fasta_bytes=before.st_size,
        identity_check='exact_stage1_sha256' if known else 'new_genomic_fingerprint; stage1_protein_search_did_not_hash_genome',
        windows=len(own),annotation_models=len(annotation))
    return dict(receipt=receipt, windows=extracted, staged=staged, annotation_features=staged_features)

def run(args):
    from .early_synteny import cache_lock
    from .stage2_resume import lock_path
    with cache_lock(lock_path(args.out)):
        return _run(args)


def _run(args):
    started=time.perf_counter();source,data=verify(args.results)
    if args.predictor!='tiberius':raise ValueError('The new stage2 command requires Tiberius')
    regions,selection,decision_rows,prediction_intervals=planned_regions(source,data,args)
    db,db_info=verify_db(args.protein_db)
    import math
    if not math.isfinite(args.identification_evalue) or args.identification_evalue<=0:raise ValueError('E-value must be positive')
    versions=tool_versions(['blastp','makeblastdb'])
    if any(not v['path'] for v in versions.values()):raise ValueError('Stage 2 requires BLAST+ blastp and makeblastdb')
    args.predictor_cache=str(Path(args.predictor_cache or Path(args.out).resolve().parent/'_genecontext_cache'/'tiberius').resolve())
    args.offline=True;args.predictor='tiberius'
    inputs=datasets(args.genomes);by_genome={r['genome_id']:r for r in inputs}
    if set(r['genome_id'] for r in regions)-set(by_genome):raise ValueError('Selected genome absent from --genomes')
    predictor_identity=validate_prediction_inputs(args,[by_genome[g] for g in sorted({r['genome_id'] for r in regions})])
    from .stage2_intervals import project_model
    from .stage2_resume import prepare as prepare_resume
    from .artifact_cache import cached
    plan_data=window_plan(source,data,args,regions,selection,prediction_intervals)
    root,run_record,input_hashes,input_stats,completed=prepare_resume(args,plan_data,by_genome,db_info,versions,predictor_identity)
    if completed is not None:
        print('STOP 2 already complete — '+str(root/'report.html'))
        return 1 if completed['errors'] else 0
    runner=TimedRunner(root);errors=[];models=[];coding={};proteins={};windows={};receipts=[];annotation_features=[]
    dump(root/'window_plan.json',plan_data)
    # Manual window/anchor overrides may appear only in later decision rows.
    columns=list(dict.fromkeys(key for row in decision_rows for key in row))
    write_tsv(root/'selection.snapshot.tsv',decision_rows,columns)
    import shutil
    decision_source=Path(selection.get('draft') or args.decisions)
    shutil.copyfile(decision_source,root/('selection.input.json' if args.draft else 'selection.input.tsv'))
    lmap={l['locus_id']:l for l in data['loci']};source_genomes={g['genome_id']:g for g in data['genomes']}
    for gid in sorted({r['genome_id'] for r in regions}):
        own=[r for r in regions if r['genome_id']==gid];entry=by_genome[gid];t0=time.perf_counter()
        receipt=dict(genome_id=gid,status='running');receipts.append(receipt)
        try:
            from .search_cache import fingerprint
            def check_inputs():
                for field,before in input_stats[gid].items():
                    if fingerprint(entry[field])!=before:raise ValueError('Input changed after hashing: '+field)
            check_inputs()
            identity=dict(schema='gcf_stage2_extraction_v1',genome_id=gid,inputs=input_hashes[gid],
                regions=own,prediction_intervals=[c for c in prediction_intervals if c['genome_id']==gid],
                source_evidence_sha256=plan_data['source_evidence_sha256'],
                implementation=run_record['identity']['implementation_sha256'])
            def build_extraction(folder):
                value=extract_genome(gid,own,entry,prediction_intervals,lmap,source_genomes,input_hashes[gid]['fasta'])
                check_inputs()
                dump(folder/'extracted.json',value)
            stored,checkpoint=cached(root/'.stage2'/'extraction',identity,build_extraction)
            value=json.loads((stored/'extracted.json').read_text())
            check_inputs()
            receipt.update(value['receipt'],checkpoint=checkpoint)
            windows.update(value['windows']);annotation_features.extend(value['annotation_features'])
            for record,cds,protein in value['staged']:
                models.append(record);coding[record['model_id']]=cds
                if protein and '*' not in protein:proteins[record['model_id']]=protein
        except (OSError,ValueError,KeyError,EOFError) as error:
            receipt.update(status='error',error=str(error));errors.append(dict(stage='extraction',genome_id=gid,error=str(error)))
            for r in own:windows.pop(r['region_id'],None)
            for c in prediction_intervals:
                if c['genome_id']==gid:windows.pop(c['region_id'],None)
        receipt['elapsed_seconds']=time.perf_counter()-t0
    write_fasta(root/'windows.fna',[(r['region_id'],windows[r['region_id']]) for r in regions if r['region_id'] in windows])
    write_fasta(root/'prediction_windows.fna',[(c['region_id'],windows[c['region_id']]) for c in prediction_intervals if c['region_id'] in windows])
    rows=[]
    for r in prediction_intervals:
        if r['region_id'] not in windows:continue
        rows.append(dict(r,locus_id=r['region_id'],source_locus_id=r['locus_id'],window_bp=len(windows[r['region_id']]),
            context_predictor_model=prediction_model(args,by_genome[r['genome_id']]),context_status='prediction_pending'))
    from .tiberius_backend import predict_batches
    (root/'context_prediction').mkdir(exist_ok=True)
    predictions,prediction_receipt=predict_batches(args,root,rows,windows,runner,errors) if rows else ({},dict(status='no_extracted_windows'))
    winseq=WindowSequence(windows);region_lookup={r['region_id']:r for r in regions}
    canonical_lookup={r['region_id']:r for r in prediction_intervals}
    for rid,values in predictions.items():
        for m in values:
            cds,protein,status,issues=coding_sequence(m,winseq)
            canonical=canonical_lookup[rid]
            for region_id in canonical['region_ids']:
                r=region_lookup[region_id];projected=project_model(m,canonical,r)
                if projected is None:continue
                record=model_record(projected,r,'tiberius',cds,protein,status,issues)
                record.update(prediction_interval_id=rid,prediction_model_id=m.id)
                models.append(record);coding[record['model_id']]=cds
                if protein and '*' not in protein:proteins[record['model_id']]=protein
    write_fasta(root/'cds.fna',sorted(coding.items()));write_fasta(root/'proteins.faa',sorted(proteins.items()))
    from .stage2_resume import protein_step
    steps={}
    step_identity=dict(schema='gcf_stage2_protein_step_v1',run=run_record['identity'],
        proteins_sha256=sha256(root/'proteins.faa'))
    identification={};id_hits=[];target_hits=[];pairs=[]
    try:
        (identification,id_hits),steps['identification']=protein_step(root,'identification',step_identity,runner,
            lambda folder:identify(runner,folder,proteins,db,db_info,args))
    except (OSError,ValueError,RuntimeError) as error:errors.append(dict(stage='identification',error=str(error)))
    try:
        queries={q['id']:q['sequence'] for q in data['queries']} if data['config']['query_type']!='dna' else {}
        target_hits,steps['target_support']=protein_step(root,'target_support',step_identity,runner,
            lambda folder:blast_pairs(runner,folder,proteins,queries,args.threads,args.identification_evalue))
        assign_groups(models)
        local,steps['neighbor_support']=protein_step(root,'neighbor_support',step_identity,runner,
            lambda folder:blast_pairs(runner,folder,proteins,proteins,args.threads,args.identification_evalue))
        pairs=neighbor_pairs(models,local,args)
    except (OSError,ValueError,RuntimeError) as error:errors.append(dict(stage='local_protein_comparison',error=str(error)))
    for m in models:
        m['identification']=identification.get(m['model_id'],dict(status='identification_error' if m['model_id'] in proteins else 'no_usable_protein',proposed_name='',near_best=[]))
        m['target_support']=sorted([h for h in target_hits if h['query']==m['model_id']],key=lambda h:(-h['bitscore'],h['subject']))
    payload=dict(schema='stage2_review_v1',version=__version__,review_id=run_record['review_id'],created=now(),
        source_review_id=data['review_id'],source_checkpoint=str(source),source_evidence_sha256=sha256(source/'evidence.json'),
        selection=selection,regions=regions,models=models,neighbor_pairs=pairs,errors=errors,annotation_features=annotation_features,
        input_genomes=[{k:g[k] for k in ['genome_id','species'] if k in g} for g in data['genomes']],
        prediction_intervals=prediction_intervals,
        query_members={q['id']:q['members'] for q in data['queries']},
        provenance=dict(genomes=receipts,prediction=prediction_receipt,protein_steps=steps,protein_db=db_info,protein_db_path=str(db),commands=runner.commands,
            input_manifest_sha256=sha256(args.genomes),
            prediction_models_sha256=sha256(args.prediction_models) if args.prediction_models else None,
            source_files_sha256={name:sha256(Path(__file__).with_name(name)) for name in ['artifact_cache.py','stage2_resume.py','tiberius_session.py','stage2.py','stage2_windows.py','annotation_context.py','stage2_intervals.py','stage2_db.py','stage2_quality.py','stage2_display.py','stage2_reconstruct.py','stage2_prior_context.py','stage2_choices.py','stage2_report.py','stage2_review.txt','stage2_neighborhood.txt','stage2_cluster.txt','stage2_review_queue.txt','tiberius_backend.py','tiberius_worker.py','sequence.py']},parameters=dict({k:getattr(args,k) for k in ['identification_identity','identification_coverage','identification_evalue','competing_ratio','threads']},window_policy=selection['window_policy'],max_window_bases=args.max_window_bases),
            elapsed_seconds=time.perf_counter()-started),
        status='completed_with_errors' if errors else 'awaiting_manual_review',
        note='No model or neighbor is automatically approved. No DNA/MSA alignment is run before STOP 2.')
    from .stage2_prior_context import collect
    payload['prior_synteny_context']=collect(payload,data,sha256(source/'evidence.json'))
    dump(root/'evidence.json',payload);dump(root/'identification_hits.json',id_hits)
    write_outputs(root,payload)
    print('STOP 2 — '+str(root/'report.html'))
    return 1 if errors else 0


def write_outputs(root,data):
    from urllib.parse import quote
    from .stage2_resume import pending_table, atomic_json
    with (root/'models.gff3').open('w') as out:
        out.write('##gff-version 3\n')
        for m in data['models']:
            prefix=f"{quote(m['genome_id'],safe='')}|{quote(m['contig'],safe='')}\t{m['source']}\t"
            out.write(f"{prefix}gene\t{m['start']+1}\t{m['end']}\t.\t{m['strand']}\t.\tID={m['model_id']};original_gene={quote(m['original_gene'],safe='')};genome_id={quote(m['genome_id'],safe='')};region_id={m['region_id']}\n")
            out.write(f"{prefix}mRNA\t{m['start']+1}\t{m['end']}\t.\t{m['strand']}\t.\tID={m['model_id']}.t1;Parent={m['model_id']}\n")
            for a,b,phase in m['cds_blocks']:
                out.write(f"{prefix}CDS\t{a+1}\t{b}\t.\t{m['strand']}\t{phase if phase is not None else '.'}\tParent={m['model_id']}.t1\n")
            for a,b in m['stop_blocks']:
                out.write(f"{prefix}stop_codon\t{a+1}\t{b}\t.\t{m['strand']}\t0\tParent={m['model_id']}.t1\n")
    annotation_dir=root/'annotations';annotation_dir.mkdir(exist_ok=True)
    native=defaultdict(list)
    for line in (root/'models.gff3').read_text().splitlines():
        if line.startswith('#'):continue
        seqid,rest=line.split('\t',1);gid,contig=seqid.split('|',1)
        from urllib.parse import unquote
        native[unquote(gid)].append(unquote(contig)+'\t'+rest)
    mapping=[]
    for gid,lines in sorted(native.items()):
        filename='annotations/'+ident(gid)+'.gff3'
        (root/filename).write_text('##gff-version 3\n'+'\n'.join(lines)+'\n')
        mapping.append(dict(genome_id=gid,gff3=filename))
    write_tsv(root/'annotations.tsv',mapping,['genome_id','gff3'])
    pending_table(root/'regions.decisions.tsv',[dict(review_id=data['review_id'],region_id=r['region_id'],decision='',notes='') for r in data['regions']],['review_id','region_id','decision','notes'])
    pending_table(root/'neighbors.decisions.tsv',[],['review_id','gene_a','gene_b','model_a','model_b','decision'])
    write_tsv(root/'regions.tsv',data['regions'],['genome_id','locus_id','region_id','contig','strand','anchor_start','anchor_end','window_start','window_end','boundary_source','window_source'])
    pending_table(root/'models.decisions.tsv',[dict(review_id=data['review_id'],model_id=m['model_id'],decision='',role='',notes='') for m in data['models']],['review_id','model_id','decision','role','notes'])
    write_tsv(root/'models.tsv',[dict(m,proposed_name=m['identification']['proposed_name'],identification_status=m['identification']['status'],issues=';'.join(m['issues'])) for m in data['models']],['model_id','region_id','genome_id','gene_group','original_model','original_gene','source','annotation_name','contig','start','end','strand','coding_status','issues','protein_length','window_edge','proposed_name','identification_status'])
    from .stage2_report import write_report
    write_report(root,data)
    atomic_json(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name not in {'seal.json','models.decisions.tsv','regions.decisions.tsv','neighbors.decisions.tsv'}})
