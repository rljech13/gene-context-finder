import copy
import unittest
import tempfile
from pathlib import Path
from genecontext.sequence import write_fasta
from genecontext.stage3 import prepare_inputs
from test_stage3 import fixture
from genecontext.stage2_choices import validate_choices
from genecontext.stage2_report import render_document

class CandidateContract(unittest.TestCase):
    def test_multiple_cds_are_comparison_candidates_not_approved_targets(self):
        data,state,*_=fixture();m=copy.deepcopy(data['models'][0])
        m.update(model_id='alternative',original_model='alternate_transcript',cds_sha256='different')
        data['models'].append(m)
        state['models']['alternative']=dict(decision='inspect',role='alignment_candidate',notes='Compare alternative CDS',target_override=False)
        self.assertEqual(validate_choices(data,state),[])
        state['models']['alternative']['role']='target'
        self.assertTrue(any('one physical gene' in x['message'].lower() for x in validate_choices(data,state)))

    def test_candidate_handoff_keeps_both_native_sequences_and_provenance(self):
        data,state,windows,cds,proteins=fixture();model=copy.deepcopy(data['models'][0]);original=model['model_id'];model.update(model_id='alternative',original_model='alternative_transcript');data['models'].append(model)
        cds['alternative']=cds[original];proteins['alternative']=proteins[original]
        state['models']['alternative']=dict(decision='inspect',role='alignment_candidate',notes='Compare alternative CDS',target_override=False)
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name,records in [('windows.fna',windows),('cds.fna',cds),('proteins.faa',proteins)]:write_fasta(root/name,records.items())
            selected=prepare_inputs(data,state,root);self.assertEqual(len(selected),3)
            candidate=next(row for row,*_ in selected if row['model_id']=='alternative')
            self.assertTrue(candidate['alignment_candidate']);self.assertFalse(candidate['target_override'])

    def test_candidates_still_require_a_reason_when_partial(self):
        data,state,*_=fixture();mid=data['models'][0]['model_id']
        state['models'][mid].update(role='alignment_candidate',notes='')
        data['models'][0]['coding_status']='partial'
        self.assertTrue(any('partial/edge' in x['message'] for x in validate_choices(data,state)))

    def test_candidate_inclusion_never_assigns_target_override(self):
        data,state,*_=fixture();mid=data['models'][0]['model_id']
        state['models'][mid].update(role='alignment_candidate',target_override=False)
        data['models'][0]['target_support']=[]
        self.assertEqual(validate_choices(data,state),[])
        self.assertFalse(state['models'][mid]['target_override'])

if __name__=='__main__':unittest.main()
