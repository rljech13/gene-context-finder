import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from genecontext.review_service import ReviewService
from genecontext.review_service_worker import write_json, verified_result

class ReviewServiceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        data=dict(review_id='review-one',stage=3,sequences=[dict(id='a',protein='MAAA',cds='ATGGCTGCTGCTTAA',group='g'),dict(id='b',protein='MAAT',cds='ATGGCTGCTACTTAA',group='g'),dict(id='c',protein='MAAA',cds='',group='other')])
        dataset=self.root/'review_data.json';write_json(dataset,data)
        self.report=dict(id='S3_fixture',review_id='review-one',stage=3,dataset=str(dataset),dataset_sha256=hashlib.sha256(dataset.read_bytes()).hexdigest(),pool_sha256='pool',paths=['/fixture/report.html'])
        self.config=self.root/'config.json';write_json(self.config,dict(state=str(self.root/'state'),reports=[self.report],runner='pinned-engine',runner_sha256='engine',engine_fingerprint='v1',max_active=2,allowed_hosts=['127.0.0.1:8794']))
        self.submissions=[];self.service=ReviewService(self.config,submit=lambda root:self.submissions.append(root))
    def test_same_selection_deduplicates_and_reload_reconnects(self):
        q={'ids':['a','b'],'group':'g'};a=self.service.align(self.report,'workspace_12345678',q);b=self.service.align(self.report,'workspace_12345678',dict(q,ids=['b','a']))
        self.assertEqual(a['job_id'],b['job_id'])
        restored=ReviewService(self.config,submit=lambda root:self.fail('Reload resubmitted a job'))
        self.assertEqual(restored.status(self.report,'workspace_12345678','g')['job']['job_id'],a['job_id'])
        self.assertEqual(restored.status(self.report,'another_workspace_123','g').get('job'),None)
    def test_source_hash_ids_groups_and_decision_checkpoint_are_enforced(self):
        for body in [{'ids':['a','a'],'group':'g'},{'ids':['a','missing'],'group':'g'},{'ids':['a','c'],'group':'g'},{'ids':['a','b'],'group':'g','decisions':{'review_id':'wrong'}}]:
            with self.assertRaises(ValueError):self.service.align(self.report,'workspace_12345678',body)
        Path(self.report['dataset']).write_text('{}')
        with self.assertRaisesRegex(ValueError,'changed'):self.service.align(self.report,'workspace_12345678',{'ids':['a','b'],'group':'g'})
    def test_session_is_bound_to_model_pool_workspace_and_report(self):
        _,_,a=self.service.session({'report':['/fixture/report.html'],'workspace':['workspace_12345678'],'pool':['pool']})
        _,_,b=self.service.session({'report':['/fixture/report.html'],'workspace':['workspace_abcdefgh'],'pool':['pool']})
        self.assertNotEqual(a,b)
        with self.assertRaises(ValueError):self.service.session({'report':['/fixture/report.html'],'workspace':['workspace_12345678'],'pool':['outdated']})
    def test_failed_attempt_is_preserved_and_explicit_retry_gets_new_job(self):
        q={'ids':['a','b'],'group':'g'};a=self.service.align(self.report,'workspace_12345678',q)
        write_json(self.root/'state/jobs'/a['job_id']/'status.json',{'status':'failed','step':'Diagnostic failure'})
        b=self.service.align(self.report,'workspace_12345678',q)
        self.assertNotEqual(a['job_id'],b['job_id']);self.assertEqual(self.service.view(a['job_id'])['status'],'failed')
    def test_result_seal_rejects_changed_alignment(self):
        root=self.root/'result';root.mkdir();write_json(root/'result.json',{'id':'example'});write_json(root/'checksums.json',{'result.json':hashlib.sha256((root/'result.json').read_bytes()).hexdigest()})
        self.assertEqual(verified_result(root)['id'],'example');(root/'result.json').write_text('{}')
        with self.assertRaises(ValueError):verified_result(root)
