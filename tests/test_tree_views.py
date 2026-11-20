import unittest
from genecontext.tree_views import rooted_tree,newick_text,views,graph
from genecontext.stage4_phylogeny import parse_newick,tip_names


class RootViews(unittest.TestCase):
    def setUp(self):
        self.tree='((A:1,B:2)98:3,(C:4,D:5)90:6,E:7);'
        self.leaves=[dict(tree_id=n,sequence_ids=[n],labels=['Species '+n],included=True) for n in 'ABCDE']

    def distances(self,text):
        adj,names=graph(parse_newick(text));out={}
        for a in adj:
            if len(adj[a])!=1:continue
            stack=[(a,None,0)]
            while stack:
                n,parent,d=stack.pop()
                if len(adj[n])==1:out[names[a],names[n]]=d
                stack.extend((k,n,d+edge[0]) for k,edge in adj[n].items() if k!=parent)
        return out

    def test_root_changes_display_without_changing_taxa_or_pair_distances(self):
        root,meta=rooted_tree(self.tree,['A','B']);rooted=newick_text(root)
        self.assertEqual(set(tip_names(root)),set('ABCDE'))
        self.assertEqual(self.distances(self.tree),self.distances(rooted))
        self.assertEqual(meta['root_edge_length'],3)
        v=views(self.tree,self.leaves,['A','B'],'Explicit family profile in this test')
        self.assertEqual(v['unrooted_newick'],self.tree)
        self.assertIn('Unrooted',v['unrooted_svg']);self.assertIn('<svg',v['rooted_svg'])
        self.assertIn('not to scale',v['rooted_svg'])
        self.assertNotIn('substitutions / site',v['rooted_svg'])
        import xml.etree.ElementTree as ET
        svg=ET.fromstring(v['unrooted_svg']);labels=[n for n in svg.iter() if 'data-tree-label' in n.attrib]
        self.assertEqual(len(labels),5)
        for anchor in ['start','end']:
            ys=sorted(float(n.attrib['y']) for n in labels if n.attrib['text-anchor']==anchor)
            self.assertTrue(all(b-a>=27.99 for a,b in zip(ys,ys[1:])))

    def test_no_forced_root_or_mixed_identical_aliases(self):
        with self.assertRaisesRegex(ValueError,'not a separate clade'):rooted_tree(self.tree,['A','C'])
        with self.assertRaises(ValueError):views(self.tree,self.leaves,['A'],'')
        self.leaves[0]['sequence_ids'].append('same_A')
        with self.assertRaisesRegex(ValueError,'Identical proteins'):views(self.tree,self.leaves,['A'],'basis')
        self.assertEqual(views(self.tree,self.leaves)['rooting'],'unrooted')

    def test_result_identity_binds_outgroup_and_its_basis(self):
        from genecontext.phylogeny_review_worker import result_identity
        def key(ids,basis):return result_identity('alignment',{'engine':'hash'},dict(outgroup_ids=ids,root_basis=basis),['code'])
        self.assertEqual(key(['A','B'],'basis'),key(['B','A'],'basis'))
        self.assertNotEqual(key(['A'],'basis'),key(['B'],'basis'))
        self.assertNotEqual(key(['A'],'basis'),key(['A'],'changed interpretation'))
        from genecontext.phylogeny_review_worker import result_identity
        self.assertNotEqual(result_identity('same-alignment',{},dict(ids=['a','b']),['code']),result_identity('same-alignment',{},dict(ids=['a','c']),['code']))


if __name__=='__main__':unittest.main()
