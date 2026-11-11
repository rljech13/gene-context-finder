"""Small offline STOP-2 workbench: models and proposed neighbor correspondences."""
import json
from pathlib import Path


def source_pairs(root, data):
    """Re-rank saved local BLAST hits using only each displayed model source.

    These are presentation evidence, not new review proposals or approvals.
    An annotation winner must never stand in for a predicted protein match.
    """
    from types import SimpleNamespace
    from .stage2 import neighbor_pairs, parse_hits
    from .sequence import read_fasta
    raw = Path(root) / 'neighbor_support' / 'blast.tsv'
    params = data.get('provenance', {}).get('parameters', {})
    if raw.exists() and raw.stat().st_size and params:
        proteins = dict(read_fasta(Path(root) / 'proteins.faa'))
        hits = parse_hits(raw, proteins, {k: len(v) for k, v in proteins.items()})
        args = SimpleNamespace(**params)
        result = {}
        for source in ['tiberius', 'annotation']:
            selected = [m for m in data['models'] if m['source'] == source]
            ids = {m['model_id'] for m in selected}
            result[source] = neighbor_pairs(selected,
                [h for h in hits if h['query'] in ids and h['subject'] in ids], args)
        return result
    lookup = {m['model_id']: m for m in data['models']}
    return {source: [p for p in data['neighbor_pairs']
                    if all(lookup[p[k]]['source'] == source for k in ['model_a', 'model_b'])]
            for source in ['tiberius', 'annotation']}


def write_report(root,data):
    from .stage2_quality import with_annotation_quality
    data=with_annotation_quality(root,data)
    from .stage2_reconstruct import load as load_reconstructions
    data=load_reconstructions(root,data)
    payload={k:data[k] for k in ['review_id','source_review_id','status','regions','models','neighbor_pairs','errors','query_members']}
    from .gene_identity import model_identities, export_identities
    payload['copy_identities']=model_identities(data)
    export_identities(root,data['models'],payload['copy_identities'],'model_id')
    target_groups={m['gene_group'] for m in data['models'] if m['overlaps_anchor']}
    payload['neighbor_pairs']=[dict(p,kind=p.get('kind') or ('target_similarity' if p['gene_a'] in target_groups or p['gene_b'] in target_groups else 'neighbor_proposal')) for p in data['neighbor_pairs']]
    payload['source_neighbor_pairs']=source_pairs(root,data)
    from .stage2_choices import constraints
    payload['model_constraints']=constraints(data)
    from .stage2_prior_context import load
    payload['prior_synteny_context']=load(root,data)
    from .stage2_display import filtered_fragments, automatic_comments, POLICY
    payload['filtered_predictions']=filtered_fragments(root,data)
    payload['automatic_comments']=automatic_comments(data,payload['filtered_predictions'])
    payload['presentation_policy']=POLICY
    payload['prediction_intervals']=data.get('prediction_intervals',[])
    payload['annotation_features']=data.get('annotation_features',[])
    payload['annotation_context_available']='annotation_features' in data
    supplement=Path(root)/'annotation_context.json'
    payload['annotation_context_external']=supplement.exists()
    if supplement.exists():
        from .sequence import sha256
        context=json.loads(supplement.read_text())
        if context.get('schema')!='stage2_annotation_context_v1' or context.get('source_evidence_sha256')!=sha256(Path(root)/'evidence.json'):
            raise ValueError('Annotation context belongs to another scientific checkpoint')
        payload['annotation_features']=context['features']
        payload['annotation_context_available']=True
    payload['prediction_summary']=dict(shared_intervals=len(data.get('prediction_intervals',[])),
        physical_genes=len({v['physical_id'] for v in payload['model_constraints'].values()}),
        model_observations=len(data['models']))
    draft_path=Path(root)/'selection.input.json'
    draft=json.loads(draft_path.read_text()) if draft_path.exists() else {}
    available={r['genome_id'] for r in data['regions']}
    all_ids=list(dict.fromkeys([g['genome_id'] for g in data.get('input_genomes',[])]+
        [p['genome_id'] for p in draft.get('pairs',[])]+list(available)+
        data.get('selection',{}).get('omitted_unsaved_genomes',[])))
    decisions={p['genome_id']:p for p in draft.get('pairs',[])}
    omitted=set(data.get('selection',{}).get('omitted_unsaved_genomes',[]))
    payload['run_context']=dict(gene=' / '.join(dict.fromkeys(r.get('gene','target') for r in data['regions'])),
        source_review_id=data['source_review_id'],genomes=[dict(genome_id=gid,
            status='processed_selected_loci' if gid in available else 'stage1_review_pending' if gid in omitted or decisions.get(gid,{}).get('committed') is False else decisions.get(gid,{}).get('action') or 'not_in_this_run') for gid in all_ids])
    template=Path(__file__).with_name('stage2_review.txt').read_text()
    template=template.replace('__NEIGHBORHOOD__',Path(__file__).with_name('stage2_neighborhood.txt').read_text())
    template=template.replace('__CLUSTER__',Path(__file__).with_name('stage2_cluster.txt').read_text())
    template=template.replace('__REVIEW_QUEUE__',Path(__file__).with_name('stage2_review_queue.txt').read_text())
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    (root/'report.html').write_text(template.replace('__DATA__',text))


def refresh(args):
    import shutil
    from .stage1 import new_output,dump
    from .sequence import sha256
    source=Path(args.results).resolve()
    for name,digest in json.loads((source/'seal.json').read_text()).items():
        path=(source/name).resolve()
        if source not in path.parents or sha256(path)!=digest:raise ValueError('Stage 2 evidence changed: '+name)
    data=json.loads((source/'evidence.json').read_text())
    if data.get('schema')!='stage2_review_v1':raise ValueError('Not a stage 2 result')
    target=Path(args.out).resolve()
    if target==source or source in target.parents or target in source.parents:raise ValueError('Use a separate presentation directory')
    root=new_output(target)
    for name in json.loads((source/'seal.json').read_text()):
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,p)
    for name in ['models.decisions.tsv','regions.decisions.tsv','neighbors.decisions.tsv']:
        if (source/name).exists():shutil.copyfile(source/name,root/name)
    if getattr(args,'genomes',None):
        from .annotation_context import collect
        context=collect(data,args.genomes);context['source_evidence_sha256']=sha256(root/'evidence.json')
        dump(root/'annotation_context.json',context)
        from .stage2_quality import collect_annotation_quality
        dump(root/'annotation_quality.json',collect_annotation_quality(data,args.genomes,sha256(source/'evidence.json')))
    write_report(root,data)
    from . import __version__
    dump(root/'presentation.json',dict(version=__version__,source=str(source),source_evidence_sha256=sha256(source/'evidence.json'),retained_review_id=data['review_id'],renderer_sha256={name:sha256(Path(__file__).with_name(name)) for name in ['stage2_report.py','stage2_quality.py','stage2_display.py','stage2_reconstruct.py','stage2_prior_context.py','stage2_choices.py','stage2_review.txt','stage2_neighborhood.txt','stage2_cluster.txt','stage2_review_queue.txt','gene_identity.py']},note='Presentation refresh only. Specific copy labels retain annotation and curated locus provenance separately from protein-database similarity. Scientific evidence, source identifiers, all-source review proposals and decisions are unchanged.'))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name not in {'seal.json','models.decisions.tsv','regions.decisions.tsv','neighbors.decisions.tsv'}})
    print('Updated STOP 2 report: '+str(root/'report.html'))
    return 0
