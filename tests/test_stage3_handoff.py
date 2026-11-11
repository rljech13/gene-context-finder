"""The final ledger must not turn unresolved or repeated observations into rejection."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.stage3_handoff import model_decisions
from tests.test_stage3 import fixture


class FinalChoiceLedger(unittest.TestCase):
    def base(self):
        data,state,*_=fixture()
        model=data['models'][0];region=data['regions'][0]
        data['models']=[model];data['regions']=[region]
        state['models']={model['model_id']:state['models'][model['model_id']]}
        state['regions']={region['region_id']:state['regions'][region['region_id']]}
        payload=dict(groups=[dict(id='group_TEST',name='TEST')],
                     models=[dict(model,group='group_TEST',species='Species a',cds_file='selected.cds',protein_file='selected.protein')])
        return data,state,payload

    def run_ledger(self,data,state,payload):
        before=json.dumps(state,sort_keys=True)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'source.evidence.json').write_text(json.dumps(data))
            (root/'decisions.snapshot.json').write_text(json.dumps(state))
            final=model_decisions(root,payload)
            with (root/'model_decision_audit.tsv').open() as handle:
                audit=list(csv.DictReader(handle,delimiter='\t'))
            with (root/'final_model_choices.tsv').open() as handle:
                selected=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual({r['model_id'] for r in selected},{m['model_id'] for m in payload['models']})
            self.assertTrue((root/'model_choice_history.tsv').is_file())
        self.assertEqual(json.dumps(state,sort_keys=True),before)
        return final,audit

    def test_automatic_size_comment_is_carried_without_replacing_review_reason(self):
        data,state,payload=self.base();rid=data['regions'][0]['region_id']
        state['automatic_comments']={rid:['Genomic span differs; review both sources.']}
        state['regions'][rid]['notes']='Reviewer decision retained'
        final,audit=self.run_ledger(data,state,payload)
        self.assertEqual(audit[0]['automatic_model_comments'],'Genomic span differs; review both sources.')
        self.assertEqual(audit[0]['saved_region_notes'],'Reviewer decision retained')
        self.assertEqual(final[0]['automatic_model_comments'],'Genomic span differs; review both sources.')

    def test_unresolved_copy_observed_in_proceed_window_remains_unresolved(self):
        data,state,payload=self.base();a=data['models'][0];region=data['regions'][0]
        other_region=dict(region,region_id='unresolved_region',anchor_start=500,anchor_end=560)
        data['regions'].append(other_region)
        q=dict(a,model_id='q_context',gene_group='physical_Q',original_gene='gene_Q',original_model='Q_transcript',
               start=500,end=560,cds_blocks=[[500,560,0]],cds_sha256='q_hash',overlaps_anchor=False)
        main=dict(q,model_id='q_main',region_id='unresolved_region',overlaps_anchor=True)
        data['models'] += [q,main]
        state['regions']['unresolved_region']=dict(decision='unresolved',saved=True,notes='Copy identity remains unresolved')
        final,audit=self.run_ledger(data,state,payload)
        unresolved=[r for r in final if not r['retained']]
        self.assertEqual(len(unresolved),1)
        self.assertEqual(unresolved[0]['status'],'Unresolved')
        self.assertEqual(unresolved[0]['alternatives'][0]['status'],'Unresolved')
        self.assertEqual(unresolved[0]['alternatives'][0]['observations'],2)
        self.assertEqual(unresolved[0]['reason'],'Copy identity remains unresolved')
        self.assertEqual(len(audit),3)

    def test_identical_alternative_is_not_retained_but_not_biologically_rejected(self):
        data,state,payload=self.base();a=data['models'][0]
        data['models'].append(dict(a,model_id='prediction',source='tiberius',original_model='g1_t1'))
        final,audit=self.run_ledger(data,state,payload)
        self.assertEqual(len(final),1)
        self.assertEqual(final[0]['retained'][0]['model_id'],a['model_id'])
        alternative=final[0]['alternatives'][0]
        self.assertEqual(alternative['status'],'Not retained')
        self.assertIn('Identical CDS retained as',alternative['reason'])
        self.assertEqual(alternative['saved_model_decision'],'')


if __name__=='__main__':unittest.main()
