import hashlib
import json
import unittest

from genecontext.synteny_reference import attach
from genecontext.structure_approval import verify


def fixture():
    handoff=dict(review_id='source',alignment_review_id='alignment',alignment_result_id='msa',tree_result_id='tree',
        decisions=dict(schema='stage2_draft_v1',review_id='source',regions={'r':dict(saved=True,decision='inspect',notes='Локус')},
            models={'a':dict(decision='inspect',role='alignment_candidate',notes='Вариант')},pairs={}))
    fingerprint=hashlib.sha256(json.dumps(dict(handoff=handoff,synteny_pairs={}),sort_keys=True,separators=(',', ':'),ensure_ascii=False).encode()).hexdigest()
    protein=hashlib.sha256(b'MAAAA').hexdigest()
    approval=dict(schema='gcf_synteny_reference_approval_v1',stage=3,source_review_id='source',alignment_review_id='alignment',
        handoff=handoff,references={'a':dict(schema='synteny_reference_approval_v1',stage=3,confirmed=True,synteny_reviewed=True,
            model_id='a',source_review_id='source',alignment_review_id='alignment',context_sha256=fingerprint,
            protein_sha256=protein,gene='TEST',genome_id='synthetic',physical_id='synthetic|locus',region_id='r')},proteins={'a':protein})
    return approval


class SyntenyReference(unittest.TestCase):
    def test_hook_is_additive_idempotent_and_requires_new_synteny_stop(self):
        legacy='<body><main>legacy choices</main></body>'
        self.assertEqual(attach(legacy),legacy)
        html='<body><main>saved DATA stays exact</main><script>// gcf-separate-synteny-stop-v2\n</script></body>'
        changed=attach(html)
        self.assertIn('<main>saved DATA stays exact</main>',changed)
        self.assertTrue(changed.index('// gcf-synteny-reference-approval-v1')>changed.index('// gcf-separate-synteny-stop-v2'))
        self.assertEqual(attach(changed),changed)
        from genecontext.interactive_review import attach as attach_review, MARKER
        self.assertIn('// gcf-synteny-reference-approval-v1',attach_review(MARKER+html,{}))

    def test_explicit_post_synteny_native_protein_is_allowed(self):
        a=fixture()
        self.assertIs(verify(a,{'a':'MAAAA'}),a)

    def test_stop2_confirmation_alone_and_missing_reference_are_insufficient(self):
        for mutate in [lambda a:a['references']['a'].update(confirmed=False),lambda a:a['references']['a'].update(synteny_reviewed=False),lambda a:a.update(references={})]:
            a=fixture();mutate(a)
            with self.assertRaises(ValueError):verify(a,{'a':'MAAAA'})

    def test_changed_decision_and_exact_protein_and_foreign_query_are_rejected(self):
        a=fixture();a['handoff']['decisions']['models']['a']['notes']='Changed'
        with self.assertRaises(ValueError):verify(a,{'a':'MAAAA'})
        with self.assertRaises(ValueError):verify(fixture(),{'a':'MGGGG'})
        with self.assertRaises(ValueError):verify(fixture(),{'query':'MAAAA'})

    def test_foreign_scope_and_unconfirmed_region_and_empty_request_are_rejected(self):
        for field in ['source_review_id','alignment_review_id']:
            a=fixture();a[field]='foreign'
            with self.assertRaises(ValueError):verify(a,{'a':'MAAAA'})
        a=fixture();a['handoff']['decisions']['regions']['r']['saved']=False
        with self.assertRaises(ValueError):verify(a,{'a':'MAAAA'})
        with self.assertRaises(ValueError):verify(fixture(),{})
