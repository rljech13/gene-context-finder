"""A bounded MAFFT + ML-tree derivative; never edits native checkpoints."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

from .review_service_worker import run as run_mafft, sha, write_json, verified_result
from .sequence import write_fasta
from .stage4_phylogeny import build_tree, tip_names, parse_newick
from .tools import Runner


def result_identity(alignment_id,settings,request,source_hashes):
    return hashlib.sha256(json.dumps([alignment_id,settings,
        sorted(request.get('ids',[])),sorted(request.get('outgroup_ids',[])),request.get('root_basis',''),source_hashes],sort_keys=True).encode()).hexdigest()[:24]


def run(job):
    job = Path(job).resolve()
    config = json.loads((job/'job.json').read_text())
    phase = config['request'].get('phase', 'combined')
    settings = config['phylogeny']
    if phase != 'alignment' and sha(settings['engine']) != settings['engine_sha256']:
        raise ValueError('Pinned IQ-TREE binary changed')
    write_json(job/'status.json', dict(status='running', step='Aligning selected candidate proteins'))
    try:
        # Reuse the existing exact protein/CDS runner; its sealed output stays intact.
        reused=None
        for candidate in (config.get('reuse_results',[]) if phase != 'alignment' else []):
            candidate=Path(candidate)
            old=json.loads((candidate/'job.json').read_text());status=json.loads((candidate/'status.json').read_text())
            if status['status']!='done' or old['dataset_sha256']!=config['dataset_sha256'] or old['runner_sha256']!=config['runner_sha256'] or old.get('phylogeny')!=settings:
                continue
            result=verified_result(status['verified_directory'])
            if result['group']!=config['request']['group'] or set(result['ids'])!=set(config['request']['ids']) or result.get('phylogeny',{}).get('status')!='computed':
                continue
            reused=(candidate,result);break
        if phase == 'tree':
            native_root = Path(config['alignment_source'])
            alignment = verified_result(native_root)
            if alignment['id'] != config['request']['alignment_result_id'] or alignment.get('phase') != 'alignment':
                raise ValueError('Tree input is not the saved alignment selection')
            ids = config['request']['ids']
            if not set(ids) <= set(alignment['ids']):
                raise ValueError('Tree sequences are outside the saved alignment')
            alignment = dict(alignment, ids=ids, protein={sid:alignment['protein'][sid] for sid in ids},
                             cds={sid:alignment['cds'][sid] for sid in ids if sid in alignment.get('cds',{})})
            mafft = native_root
        else:
            mafft = (reused[0] if reused else job)/'mafft'
            if not reused:
                mafft.mkdir()
                alignment_config={k:v for k,v in config.items() if k not in ['phylogeny','reuse_results']}
                alignment_config['request']={k:v for k,v in config['request'].items() if k in ['ids','group','decisions']}
                write_json(mafft/'job.json', alignment_config)
                run_mafft(mafft)
            native_root = Path(json.loads((mafft/'status.json').read_text())['verified_directory'])
            alignment = verified_result(native_root)
        pool = json.loads(Path(config['dataset']).read_text())
        by_id = {s['id']:s for s in pool['sequences']}
        if alignment['review_id'] != pool['review_id'] or set(alignment['ids']) != set(config['request']['ids']):
            raise ValueError('Alignment does not match requested review/model set')
        for sid in alignment['ids']:
            if alignment['protein'][sid].replace('-', '') != by_id[sid]['protein']:
                raise ValueError('Aligned sequence changed: '+sid)
        target = job/'derived';target.mkdir()
        folder = target/'phylogeny';runner=Runner(target)
        if phase == 'alignment':
            tree = dict(status='not_requested', reason='Build the tree after reviewing this saved alignment')
        else:
            write_json(job/'status.json', dict(status='running', step='IQ-TREE: ModelFinder and 1000 UFBoot replicates'))
            if reused and reused[1]['protein'] != alignment['protein']:
                reused = None
            if reused:
                folder.mkdir()
                tree=dict(reused[1]['phylogeny']);tree['reused_ml_result_sha256']=sha(reused[0]/'result.json')
                tree['source_commands']=tree.get('commands',[])
                unrooted=tree.get('unrooted_newick') or tree['newick']
                tree['source_ml_elapsed_seconds']=tree.pop('elapsed_seconds',None);tree['ml_reused']=True
                for key in ['svg','newick','unrooted_newick','rooted_newick','unrooted_svg','rooted_svg',
                            'unrooted_phylogram_svg','rooted_phylogram_svg','rooting','root_basis',
                            'outgroup_sequence_ids','outgroup_tree_ids','root_metadata']:
                    tree.pop(key,None)
                for name in ['tree.aligned.faa']:
                    shutil.copyfile(Path(json.loads((reused[0]/'status.json').read_text())['verified_directory'])/'phylogeny'/name,folder/name)
                (folder/'ml.treefile').write_text(unrooted)
            else:
                tree = build_tree(folder, alignment['protein'], {sid:by_id[sid]['label'] for sid in alignment['ids']},
                    runner, threads=2, bootstrap=1000, occupancy=.5, engine=settings['engine'])
            if tree['status'] == 'computed':
                newick = (folder/'ml.treefile').read_text()
                expected = {leaf['tree_id'] for leaf in tree['leaves'] if leaf['included']}
                if set(tip_names(parse_newick(newick))) != expected:
                    raise ValueError('Computed tree leaf set differs from selected alignment')
                from .tree_views import views
                displays=views(newick,tree['leaves'],config['request'].get('outgroup_ids',[]),config['request'].get('root_basis',''))
                tree.update(displays,svg=displays.get('rooted_svg',displays['unrooted_svg']),newick=newick)
        tree['commands'] = runner.commands
        tree['sequence_set_sha256'] = hashlib.sha256(json.dumps(sorted(
            (sid,hashlib.sha256(by_id[sid]['protein'].encode()).hexdigest()) for sid in alignment['ids'])).encode()).hexdigest()
        result = dict(alignment, phylogeny=tree, workflow='gcf_alignment_first_v1', phase=phase)
        # Alignment id is retained and a distinct result digest binds this derivative.
        result['alignment_id'] = alignment['id']
        from . import tree_views
        result['id'] = result_identity(alignment['id'],dict(settings,phase=phase),config['request'],[sha(__file__),sha(tree_views.__file__),config['request'].get('notes','')])
        write_json(target/'result.json', result)
        write_fasta(target/'protein.aligned.faa', result['protein'].items())
        write_fasta(target/'cds.aligned.fna', result.get('cds',{}).items())
        for name in ['proteins.native.faa','cds.native.fna','selected_models.json','source_history.json']:
            shutil.copyfile(native_root/name,target/name)
        if phase == 'alignment':
            shutil.copyfile(job/'model_selection.json', target/'model_selection.json')
        elif phase == 'tree':
            shutil.copyfile(native_root/'model_selection.json', target/'model_selection.json')
            for source, name in [('protein.aligned.faa','source_alignment.protein.faa'),('cds.aligned.fna','source_alignment.cds.fna')]:
                shutil.copyfile(native_root/source, target/name)
            write_fasta(target/'proteins.native.faa', ((sid,by_id[sid]['protein']) for sid in result['ids']))
            write_fasta(target/'cds.native.fna', ((sid,by_id[sid]['cds']) for sid in result['ids'] if by_id[sid].get('cds')))
            selected_metadata = json.loads((native_root/'selected_models.json').read_text())
            write_json(target/'selected_models.json', {sid:selected_metadata[sid] for sid in result['ids']})
            write_json(target/'tree_selection.json', dict(alignment_result_id=alignment['id'], ids=result['ids'],
                outgroup_ids=config['request'].get('outgroup_ids',[]), model_selection_unchanged=True))
        write_json(target/'tree.json',tree)
        if tree.get('svg'):(target/'tree.svg').write_text(tree['svg'])
        if tree.get('newick'):(target/'tree.nwk').write_text(tree['newick'])
        write_json(target/'provenance.json', dict(dataset_sha256=config['dataset_sha256'],
            source_review_id=pool.get('source_review_id'), review_id=pool['review_id'],
            iqtree_sha256=settings['engine_sha256'], mafft_result_sha256=sha(native_root/'result.json'), comparison_phase=phase,
            selected_ids=result['ids'], biological_approval_asserted=False))
        # Only compact reusable exports are released, not engine working trees.
        checks={p.name:sha(p) for p in target.iterdir() if p.is_file()}
        write_json(target/'checksums.json',checks)
        verified_result(target)
        with zipfile.ZipFile(job/'alignment.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
            for p in target.iterdir():
                if p.is_file():archive.write(p,p.name)
        write_json(job/'result.json',result)
        write_json(job/'status.json',dict(status='done',step='Saved selection and alignment verified' if phase == 'alignment' else 'Alignment and preliminary tree verified',
            result_id=result['id'],result_sha256=sha(job/'result.json'),
            archive_sha256=sha(job/'alignment.zip'),verified_directory=str(target)))
    except Exception as error:
        write_json(job/'status.json',dict(status='failed',step=str(error)))
        raise
