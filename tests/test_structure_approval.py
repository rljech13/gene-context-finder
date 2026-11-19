import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.structure_approval import verify,issue


def approved_sequences(sequences):
    """Explicit synthetic STOP3 approvals for offline protocol fixtures only."""
    return dict(schema='gcf_final_reference_approval_v1',stage=3,
        source_seal_sha256='0'*64,spec=dict(review_id='synthetic',set_id='synthetic',
            buckets=[dict(key=sid,model_ids=[sid]) for sid in sequences]),
        decisions=dict(schema='stage3_final_choices_v1',review_id='synthetic',set_id='synthetic',
            choices={sid:dict(primary=sid,reason='Explicit synthetic final reference approval',reviewed=True) for sid in sequences}),
        proteins={sid:hashlib.sha256(seq.encode()).hexdigest() for sid,seq in sequences.items()})


class FinalReferenceGate(unittest.TestCase):
    def test_no_saved_candidate_or_alignment_can_authorize_export(self):
        from genecontext.af3_requests import prepare
        from genecontext.af3_campaign import prepare as campaign
        for producer in [prepare,campaign]:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)/'out'
                for approval in [None,dict(schema='stage2_draft_v1'),dict(schema='saved_model_selection_v1')]:
                    with self.assertRaisesRegex(ValueError,'final reference-model approval'):producer({'m':'MAAA'},root,approval=approval)
                    self.assertFalse(root.exists())

    def test_unreviewed_nonprimary_and_changed_proteins_are_rejected(self):
        seqs={'m':'MAAA'};approval=approved_sequences(seqs)
        self.assertEqual(verify(approval,seqs)['stage'],3)
        with self.assertRaisesRegex(ValueError,'exact protein'):verify(approval,{'m':'MGGG'})
        with self.assertRaisesRegex(ValueError,'Only explicitly'):verify(approval,{'other':'MAAA'})
        approval['decisions']['choices']['m']['reviewed']=False
        with self.assertRaisesRegex(ValueError,'Only explicitly'):verify(approval,seqs)

    def test_native_stop3_approval_exports_only_reviewed_primary(self):
        from test_stage4 import fixture
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();fixture(source,True)
            approval=issue(source,root/'approval.json')
            self.assertTrue(approval['proteins'])
            allowed={v['primary'] for v in approval['decisions']['choices'].values() if v['reviewed']}
            self.assertEqual(set(approval['proteins']),allowed)
            with self.assertRaisesRegex(ValueError,'new final approval'):issue(source,root/'approval.json')

    def test_preview_does_not_export_unapproved_server_jobs(self):
        from test_stage4 import fixture
        from genecontext.cli import parser
        from genecontext.stage4 import run
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();fixture(source,False)
            args=parser().parse_args(['stage4','--results',str(source),'--out',str(root/'out'),
                '--tree','none','--preview','--structure-provider','alphafold-server'])
            with self.assertRaisesRegex(ValueError,'Only explicitly'):run(args)
            self.assertFalse((root/'out').exists())


if __name__=='__main__':unittest.main()
