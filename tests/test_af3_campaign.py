import json
from pathlib import Path
import tempfile
import unittest

from genecontext.af3_campaign import prepare as prepare_with_approval
from test_structure_approval import approved_sequences

def prepare(sequences,*args,**kwargs):
    return prepare_with_approval(sequences,*args,approval=approved_sequences(sequences),**kwargs)


class Campaign(unittest.TestCase):
    def test_every_unique_protein_is_queued_once_in_balanced_bounded_batches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);alphabet='ACDEFGHIKLMNPQRSTVWY'
            seqs={str(i):'MAAA'+alphabet[i//20]+alphabet[i%20] for i in range(65)}
            seqs['same_protein_other_species']=seqs['0']
            meta={sid:dict(gene='A' if sid.isdigit() and int(sid)<33 else 'B',species=sid) for sid in seqs}
            plan=prepare(seqs,root/'campaign',meta,first_batch_size=5,batch_size=30)
            self.assertEqual(plan['counts'],dict(selected_sequences=66,unique_requested=65,available=0,awaiting_submission=65,unsupported=0))
            self.assertEqual([u['jobs'] for u in plan['uploads']],[5,30,30])
            names=[]
            for upload in plan['uploads']:
                jobs=json.loads((root/'campaign'/upload['file']).read_text());names.extend(j['name'] for j in jobs)
            self.assertEqual(len(set(names)),65)
            first=json.loads((root/'campaign'/plan['uploads'][0]['file']).read_text())
            self.assertTrue(first[0]['name'].startswith('GCF_A_'));self.assertTrue(first[1]['name'].startswith('GCF_B_'))
            alias=next(j for j in plan['jobs'] if '0' in j['source_sequence_ids'])
            self.assertEqual(alias['source_sequence_ids'],['0','same_protein_other_species'])
            self.assertEqual(plan['display_policy'],'one_best_structure_per_protein')

    def test_unsupported_sequences_are_explicit_and_never_sent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);plan=prepare({'good':'MAAA','unknown':'MAXA','oversize':'M'*5001},root/'campaign')
            self.assertEqual(plan['counts']['awaiting_submission'],1)
            self.assertEqual(plan['counts']['unsupported'],2)

    def test_completed_full_sequence_is_reused_and_changed_results_rejected(self):
        try:import gemmi
        except ImportError:self.skipTest('Optional gemmi not installed')
        from test_alphafold_server import fixture,archive
        from genecontext.alphafold_server import import_result
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);_,files=fixture(root)
            import_result(root/'requests/requests.json',archive(root,files),root/'result')
            plan=prepare({'new_alias':'MAAAA','new_protein':'MGGGG'},root/'campaign',completed=[root/'result'])
            self.assertEqual(plan['counts']['available'],1)
            self.assertEqual(plan['counts']['awaiting_submission'],1)
            jobs=json.loads((root/'campaign'/plan['uploads'][0]['file']).read_text())
            self.assertEqual(jobs[0]['sequences'][0]['proteinChain']['sequence'],'MGGGG')
            (root/'result/model.pdb').write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):prepare({'alias':'MAAAA'},root/'bad',completed=[root/'result'])
            self.assertFalse((root/'bad').exists())


if __name__=='__main__':unittest.main()
