from copy import deepcopy
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.alignment_review import build, stage1_pool, validate_selection


class AlignmentFirst(unittest.TestCase):
    def test_candidate_review_does_not_inherit_approval_and_binds_new_namespace(self):
        with TemporaryDirectory() as tmp:
            original={'sequences':[dict(id='m',protein='MAL',cds='ATGGCTCTG',group='g',reference=False,
                label='species · model',metadata={'source':'annotation'})], 'selected':['m'], 'visible':['m']}
            history=dict(kind='legacy_stop2_reuse',source_review_id='old',copied_approvals=False)
            with patch('genecontext.alignment_review.legacy_pool',return_value=(deepcopy(original),history)):
                root,pool=build('source',Path(tmp)/'out',legacy=True)
            self.assertNotEqual(pool['review_id'],'old')
            self.assertEqual(pool['sequences'],original['sequences'])
            self.assertFalse(json.loads((root/'manifest.json').read_text())['biological_approval_asserted'])
            self.assertFalse(json.loads((root/'manifest.json').read_text())['tiberius_run'])
            self.assertIn('alignment and gene tree',(root/'report.html').read_text())
            choice=dict(schema='alignment_model_selection_v1',review_id=pool['review_id'],source_review_id='old',
                ids=['m'],biological_approval_asserted=False)
            self.assertEqual(validate_selection(pool,choice)[0]['id'],'m')
            for key,value in [('review_id','old'),('source_review_id','other'),('ids',['unknown']),
                              ('ids',['m','m']),('biological_approval_asserted',True)]:
                with self.assertRaises(ValueError):validate_selection(pool,dict(choice,**{key:value}))

    def test_stop1_exact_products_and_queries_without_prediction_or_genomic_fasta(self):
        with TemporaryDirectory() as tmp:
            source=Path(tmp);folder=source/'raw/species.annotation';folder.mkdir(parents=True)
            seq='MADE';mapping=dict(gene_id='gene_G',model_id='rna_G1',sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
            (folder/'products.json').write_text(json.dumps({'products':[{'search_id':'P1','mappings':[mapping]}]}))
            (folder/'matched_products.faa').write_text('>P1\n'+seq+'\n')
            for name in ['evidence.json','seal.json']:(source/name).write_text('{}')
            data=dict(review_id='r1',loci=[dict(locus_id='l',annotation_gene_id='gene_G',discovery_engine='blastp')],
                config=dict(gene='G'),queries=[dict(id='q',sequence=seq,members=[dict(species='reference')])])
            regions=[dict(locus_id='l',genome_id='species',gene='G',notes='Original manual comment')]
            with patch('genecontext.alignment_review.verify',return_value=(source,data)),patch(
                'genecontext.stage2.selected_regions',return_value=(regions,{'decision_source_sha256':'exact'},[])):
                pool,provenance=stage1_pool(source,SimpleNamespace())
            self.assertEqual([s['protein'] for s in pool['sequences']],[seq,seq])
            self.assertTrue(all(s['cds']=='' for s in pool['sequences']))
            self.assertEqual(provenance['tiberius'],'not_run')
            self.assertEqual(pool['sequences'][0]['metadata']['stop1_comment'],'Original manual comment')
            # Even with an otherwise valid mapping, wrong saved protein identity fails closed.
            mapping['sequence_sha256']='wrong'
            (folder/'products.json').write_text(json.dumps({'products':[{'search_id':'P1','mappings':[mapping]}]}))
            with patch('genecontext.alignment_review.verify',return_value=(source,data)),patch(
                'genecontext.stage2.selected_regions',return_value=(regions,{},[])),self.assertRaisesRegex(ValueError,'identity differs'):
                stage1_pool(source,SimpleNamespace())


if __name__=='__main__':unittest.main()
