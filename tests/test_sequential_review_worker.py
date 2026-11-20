import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from genecontext.phylogeny_review_worker import run
from genecontext.review_service_worker import sha, write_json, verified_result


class SequentialWorker(unittest.TestCase):
    def test_alignment_does_not_build_tree_and_tree_does_not_realign(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);dataset=root/'pool.json';engine=root/'iqtree';engine.write_bytes(b'pinned test engine')
            write_json(dataset,dict(review_id='r',source_review_id='source',sequences=[dict(id='a',label='A',protein='MA',cds='',group='g'),dict(id='b',label='B',protein='MT',cds='',group='g'),dict(id='c',label='C',protein='MG',cds='',group='g')]))
            cfg=dict(dataset=str(dataset),dataset_sha256=sha(dataset),runner='test',runner_sha256='test',review_id='r',report='test',phylogeny=dict(engine=str(engine),engine_sha256=sha(engine)),request=dict(phase='alignment',ids=['a','b','c'],group='g',notes='exact comment'))
            first=root/'first';first.mkdir();write_json(first/'job.json',cfg);write_json(first/'model_selection.json',dict(ids=['a','b','c'],notes='exact comment',biological_approval_asserted=False))
            def mafft(job):
                native=job/'native';native.mkdir();r=dict(id='mafft-id',review_id='r',ids=['a','b','c'],group='g',protein=dict(a='MA',b='MT',c='MG'),cds={},mafft_version='fixture')
                write_json(native/'result.json',r)
                for f in ['proteins.native.faa','cds.native.fna','selected_models.json','source_history.json']:(native/f).write_text('preserved source '+f)
                write_json(native/'selected_models.json',{sid:{'id':sid} for sid in ['a','b','c']});write_json(native/'checksums.json',{p.name:sha(p) for p in native.iterdir()});write_json(job/'status.json',dict(verified_directory=str(native)));write_json(job/'result.json',r)
            with patch('genecontext.phylogeny_review_worker.run_mafft',side_effect=mafft) as align,patch('genecontext.phylogeny_review_worker.build_tree',side_effect=AssertionError('Tree must not run')):
                run(first);self.assertEqual(align.call_count,1)
            result=verified_result(first/'derived');self.assertEqual(result['phase'],'alignment');self.assertEqual(result['phylogeny']['status'],'not_requested')
            second=root/'second';second.mkdir();cfg.update(alignment_source=str(first/'derived'),request=dict(phase='tree',ids=['a','b'],group='g',alignment_result_id=result['id']))
            write_json(second/'job.json',cfg)
            with patch('genecontext.phylogeny_review_worker.run_mafft',side_effect=AssertionError('Must use saved alignment')),patch('genecontext.phylogeny_review_worker.build_tree',return_value=dict(status='not_computed',reason='Two-tip fixture')) as tree:
                run(second);self.assertEqual(tree.call_count,1)
            final=verified_result(second/'derived');self.assertEqual(final['alignment_id'],result['id']);self.assertEqual(final['protein'],{k:result['protein'][k] for k in ['a','b']});self.assertEqual(json.loads((second/'derived/selected_models.json').read_text()).keys(),{'a','b'});self.assertEqual((second/'derived/source_alignment.protein.faa').read_bytes(),(first/'derived/protein.aligned.faa').read_bytes());self.assertEqual((second/'derived/model_selection.json').read_bytes(),(first/'model_selection.json').read_bytes())
