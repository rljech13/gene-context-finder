import json
import unittest
from genecontext.legacy_review import encoded, backtranslate, import_analysis

class CompletedReviewImport(unittest.TestCase):
    def test_html_safe_json_preserves_comments(self):
        value={'comment':'<gene> & </script> Причина выбора'}
        text=encoded(value)
        self.assertNotIn('<',text)
        self.assertEqual(json.loads(text),value)

    def test_saved_alignment_uses_exact_native_codons(self):
        member={'id':'m1','protein':'MA','cds':'ATGGCTTAA'}
        self.assertEqual(backtranslate('-M-A',member),'---ATG---GCT')
        with self.assertRaises(ValueError):backtranslate('MC',member)

    def test_missing_tree_is_not_scientific_inference(self):
        members={'m1':{'id':'m1','protein':'MA','cds':'ATGGCT'},'m2':{'id':'m2','protein':'MA','cds':'ATGGCC'}}
        result=import_analysis('gene',{'protein':{'m1':'M-A','m2':'M-A'}},members,'review',{})
        self.assertEqual(result['alignment']['protein'],{'m1':'M-A','m2':'M-A'})
        self.assertEqual(result['alignment']['cds']['m2'],'ATG---GCC')
        tree=result['tree']['phylogeny']
        self.assertEqual(tree['status'],'not_computed')
        self.assertIsNone(tree['rooted_newick'])
        self.assertEqual(tree['leaves'][0]['sequence_ids'],['m1','m2'])
        self.assertFalse(result['alignment']['provenance']['new_alignment_calculation'])

    def test_mismatching_saved_member_rejected(self):
        with self.assertRaises(ValueError):import_analysis('gene',{'protein':{'m1':'MA','m2':'MC'}},{'m1':{'protein':'MA'},'m2':{'protein':'MA'}},'review',{})

if __name__=='__main__':unittest.main()
