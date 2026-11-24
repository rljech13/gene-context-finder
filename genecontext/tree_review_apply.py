"""Apply a validated Stage4 tree review to a NEW sealed dossier; no searches."""
import csv
import hashlib
import json
from pathlib import Path
import shutil

from .sequence import sha256
from .stage1 import new_output, dump, now
from .tree_placement import METHOD, Tree, attach_stage4


def copy_sealed(source, target):
    source = Path(source).resolve(); seal = json.loads((source/'seal.json').read_text())
    for name, expected in seal.items():
        path = (source/name).resolve()
        if not path.is_relative_to(source) or not path.is_file() or sha256(path) != expected:
            raise ValueError('Stage4 source seal differs: '+name)
    target = new_output(Path(target).resolve())
    for name in seal:
        dest = target/name; dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, dest)
    shutil.copyfile(source/'seal.json', target/'source.stage4.seal.json')
    return target


def publish(root, data, source, operation):
    for filename, role in [('final_table.tsv','primary_target'),('special_cases.tsv','special_case')]:
        old = root/filename
        rows = []
        if old.exists():
            with old.open(newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t'); fields = list(reader.fieldnames or []); rows = list(reader)
            lookup = {r['model_id']:r for r in data['rows'] if r['classification']==role}
            extra = ['tree_score','tree_proposed_score','tree_reviewed','tree_assessment_status','tree_quartet_coverage']
            for name in extra:
                if name not in fields: fields.append(name)
            for row in rows:
                model = lookup.get(row['model_id'])
                if model:
                    for name in extra: row[name] = model.get(name) if model.get(name) is not None else 'NA'
            with old.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t');writer.writeheader();writer.writerows(rows)
    dump(root/'evidence.json', data)
    text = json.dumps(data,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    (root/'report.html').write_text(Path(__file__).with_name('stage4_report.txt').read_text().replace('__DATA__',text))
    (root/'START_HERE.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Stage4</title><p><a href="report.html">Stage4 dossier</a> · <a href="TREE_REVIEW.html">Tree review</a></p><p>Run python3 START_REVIEW.py for interactive tree/outgroup review.</p>')
    old = json.loads((root/'manifest.json').read_text())
    old.update(created=now(), tree_review_operation=operation, tree_placement_method=METHOD,
               source_stage4_seal_sha256=sha256(Path(source)/'seal.json'), biological_approval_asserted=False,
               tree_review_code_sha256={name:sha256(Path(__file__).with_name(name)) for name in ['tree_placement.py','tree_review_apply.py','tree_review.py','tree_review_service.py']})
    dump(root/'manifest.json', old)
    dump(root/'seal.json', {str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json' and not any(x.startswith('.') for x in p.relative_to(root).parts)})


def prepare(args):
    source = Path(args.results).resolve(); data = json.loads((source/'evidence.json').read_text())
    if data.get('schema') != 'stage4_dossier_v1':
        raise ValueError('Expected a native Stage4 dossier; study-specific adapters use the same core separately')
    root = copy_sealed(source, args.out)
    attach_stage4(root, data['trees'], data['rows'], args)
    publish(root, data, source, 'prepare_automatic_tree_review')
    return root


def apply(args):
    from .tree_review_service import prepare_overlay
    source = Path(args.results).resolve()
    inputs = json.loads((source/'tree_inputs.json').read_text())
    returned = json.loads(Path(args.reviews).read_text())
    overlay = prepare_overlay(returned, inputs)
    data = json.loads((source/'evidence.json').read_text())
    if data.get('schema') != 'stage4_dossier_v1':
        raise ValueError('Expected a native Stage4 dossier')
    root = copy_sealed(source, args.out)
    by_group = {t['group']:t for t in data['trees']}
    for panel, review in overlay['panels'].items():
        for gid, accepted in review['groups'].items():
            tree = by_group[gid]; folder = root/tree['folder']
            newick = accepted['newick']; display = accepted.get('display_newick') or newick
            tips = set(Tree(newick).tips)
            for leaf in tree['leaves']:leaf['included'] = leaf['tree_id'] in tips
            tree.update(status='computed', unrooted_newick=newick, rooting='outgroup' if accepted.get('outgroup_ids') else 'unrooted',
                        root_basis=accepted.get('root_basis',''), outgroup_tree_ids=accepted.get('outgroup_tree_ids',[]),
                        tree_reviewed=bool(accepted['approved']), tree_review_signature=accepted['signature'])
            tree['mask']=dict(original_columns=len(next(iter(inputs['panels'][panel]['groups'][gid]['rows'].values()))),
                retained_columns=len(accepted['columns']),retained_columns_1based=[i+1 for i in accepted['columns']],
                review_mask=accepted['mask'],excluded_all_gap_taxa=[])
            if accepted.get('outgroup_ids'):tree['rooted_newick']=display
            else:tree.pop('rooted_newick',None)
            (folder/'ml.treefile').write_text(display+'\n');(folder/'unrooted.nwk').write_text(newick+'\n')
            from .tree_views import unrooted_svg
            from .stage4_phylogeny import tree_svg
            (folder/'unrooted.svg').write_text(unrooted_svg(newick,tree['leaves'],True))
            (folder/'tree.svg').write_text(tree_svg(display,tree['leaves']) if accepted.get('outgroup_ids') else unrooted_svg(newick,tree['leaves'],True))
            if accepted.get('alignment_fasta'):
                (folder/'reviewed.tree.aligned.faa').write_text(accepted['alignment_fasta'])
                tree['used_alignment']='reviewed.tree.aligned.faa'
            elif accepted.get('sourceId')=='rebuilt':
                tree.pop('used_alignment',None)
            placement = accepted['placement'];dump(folder/'placement.json',placement)
            tree['placement']=dict(method=METHOD,file=str((folder/'placement.json').relative_to(root)),
                                   reference_sha256=placement['reference_sha256'],tree_sha256=placement['tree_sha256'])
            for row in data['rows']:
                if row['group']!=gid:continue
                result=review['models'][row['model_id']];raw=result['placement']
                if raw.get('protein_sha256')!=row['protein_sha256']:
                    raise ValueError('Reviewed model protein differs')
                row.update(tree_placement=raw,tree_score=result['points'],tree_proposed_score=result['proposed_points'],
                           tree_reviewed=bool(accepted['approved']),tree_review_signature=accepted['signature'],
                           tree_assessment_status=raw['status'] if accepted['approved'] else 'awaiting_tree_review',
                           tree_assessed_quartets=raw['counts'].get('assessed',0),tree_quartet_coverage=raw.get('coverage'))
    dump(root/'tree_review.applied.json', overlay)
    dump(root/'tree_initial_reviews.json', returned)
    publish(root, data, source, 'apply_signature_bound_tree_review')
    return root
