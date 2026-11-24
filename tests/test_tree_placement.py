"""Automatic placement must not depend on chosen roots, IDs, or duplicate isoforms."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.tree_placement import assess, read_lineages, triples, Tree, lineage_newick


class Placement(unittest.TestCase):
    reference = '((A,B),(C,D));'
    tree = '((a,b)80,(c,d)80);'

    def models(self):
        return [dict(model_id='model_'+taxon, tree_id=taxon.lower(), taxon=taxon, primary=True) for taxon in 'ABCD']

    def result(self, tree=None, models=None, reference=None):
        return assess(tree or self.tree, models or self.models(), species_tree=reference or self.reference, details=True)

    def test_known_concordance_and_conflict(self):
        result=self.result()
        self.assertEqual([r['points'] for r in result['models'].values()],[12.]*4)
        bad=self.result('((a,c)100,(b,d)100);')
        self.assertEqual([r['points'] for r in bad['models'].values()],[0.]*4)
        self.assertTrue(all(r['counts']['discordant']==1 for r in bad['models'].values()))

    def test_root_and_sibling_order_do_not_change_score(self):
        original=self.result()
        rerooted=self.result('(a,(b,(d,c)80));')
        for mid in original['models']:
            self.assertEqual(original['models'][mid]['points'],rerooted['models'][mid]['points'])
            self.assertEqual(original['models'][mid]['counts'],rerooted['models'][mid]['counts'])

    def test_polytomy_and_missing_support_are_na_not_zero(self):
        for tree in ['(a,b,c,d);','((a,b),(c,d));','((a,b)nan,(c,d));','((a,b)80:0,(c,d)80:0);']:
            result=self.result(tree)
            self.assertTrue(all(r['points'] is None for r in result['models'].values()),tree)
        self.assertTrue(all(r['points'] is None for r in self.result(reference='(A,B,C,D);')['models'].values()))

    def test_support_fields_never_guessed(self):
        multi='((a,b)70/99,(c,d)70/99);'
        self.assertIsNone(self.result(multi)['models']['model_A']['points'])
        r=assess(multi,self.models(),species_tree=self.reference,support_field=1)
        self.assertAlmostEqual(r['models']['model_A']['points'],14.85)

    def test_model_and_gene_names_are_not_scoring_inputs(self):
        models=self.models(); renamed=[]
        for i,m in enumerate(reversed(models)):
            renamed.append(dict(m,model_id='unseen_gene_isoform_'+str(i)))
        by_taxon={r['taxon']:r['points'] for r in self.result(models=renamed)['models'].values()}
        self.assertEqual(by_taxon,{t:12. for t in 'ABCD'})

    def test_same_species_isoforms_do_not_become_independent_taxa(self):
        models=self.models();models.append(dict(models[0],model_id='alt_A',primary=False))
        r=self.result(models=models)
        self.assertEqual(r['models']['alt_A']['points'],12.)
        for mid in ['model_B','model_C','model_D']:
            self.assertEqual(r['models'][mid]['counts']['assessed'],1)
        # An unclassified second copy cannot silently be picked as the best peer.
        models[-1].update(tree_id='x',primary=True)
        r=self.result('((a,b)80,(x,(c,d)80));',models)
        self.assertIsNone(r['models']['model_B']['points'])
        self.assertEqual(r['models']['model_A']['points'],12.)

    def test_identical_multi_species_leaf_is_unresolved(self):
        models=self.models();models[1]['tree_id']='a'
        r=self.result('(a,c,d);',models)
        self.assertTrue(all(m['points'] is None for m in r['models'].values()))

    def test_excluded_focal_and_unknown_taxon_have_explicit_reasons(self):
        models=self.models();models.append(dict(model_id='absent',tree_id='nope',taxon='A'))
        models.append(dict(model_id='unknown',tree_id='a',taxon='Unknown organism',primary=False))
        r=self.result(models=models)
        self.assertEqual(r['models']['absent']['reason'],'model_not_in_tree')
        self.assertEqual(r['models']['unknown']['reason'],'species_not_in_reference_scaffold')

    def test_alias_conflicts_duplicate_tips_and_negative_lengths_rejected(self):
        models=self.models();models[0]['protein_sha256']='one';models.append(dict(models[0],model_id='fake',protein_sha256='two'))
        with self.assertRaises(ValueError):self.result(models=models)
        for tree in ['((a,a)80,(c,d));','((a,b)80:-1,(c,d));','((a,b)80,(c,d))']:
            with self.assertRaises(ValueError):self.result(tree)

    def test_reference_taxonomy_preserves_unresolved_relationships(self):
        lineages={'A':('root','first'),'B':('root','first'),'C':('root','other'),'D':('root','other')}
        r=assess(self.tree,self.models(),lineages=lineages)
        self.assertEqual(r['models']['model_A']['points'],12.)
        self.assertEqual(r['reference_kind'],'taxonomic_scaffold')
        self.assertEqual(len(read_lineages()),18)
        quoted=Tree("(('Species a','Species b'),(C,D));",normalize=True)
        self.assertIn('Species a',quoted.bits)

    def test_large_panel_sampling_is_bounded_and_deterministic(self):
        taxa=['S'+str(n) for n in range(70)]
        one,total=triples(taxa,'focal');two,_=triples(taxa,'focal')
        self.assertEqual(one,two);self.assertEqual(len(one),2000)
        self.assertEqual(len(set(one)),2000);self.assertGreater(total,len(one))


if __name__=='__main__':unittest.main()
