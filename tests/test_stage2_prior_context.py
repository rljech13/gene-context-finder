import unittest
from genecontext.stage2_prior_context import collect

class PriorContext(unittest.TestCase):
    def test_retains_wider_evidence_separate_from_prediction_window(self):
        data=dict(source_review_id='s1',source_evidence_sha256='hash',regions=[dict(region_id='r',genome_id='g',locus_id='l',window_start=100,window_end=200)])
        wider=dict(neighbors=[dict(name='NPDC1',start=0,end=50)],comparisons=[dict(shared_anchors=1)])
        prior=dict(review_id='s1',loci=[dict(genome_id='g',locus_id='l',early_synteny=wider)])
        result=collect(data,prior,'hash')
        self.assertEqual(result['regions']['r'],wider)
        self.assertEqual(data['regions'][0]['window_start'],100)
        with self.assertRaisesRegex(ValueError,'fingerprint'):collect(data,prior,'other')
        prior['review_id']='other'
        with self.assertRaisesRegex(ValueError,'another STOP-1'):collect(data,prior,'hash')
