"""Naming guards for similar proteins at distinct, independently named loci."""
from copy import deepcopy
import unittest

from genecontext.gene_identity import identity, model_identities


class GeneIdentity(unittest.TestCase):
    def chicken(self, **changes):
        result = dict(model_id='chicken_like', genome_id='Gallus_gallus_GCF_016699485.2',
            source='annotation', annotation_name='LOC101748883', original_gene='gene-LOC101748883',
            contig='NC_052532.1', strand='+', cds_blocks=[[75862830,75862914,0],
            [75870577,75870700,0],[75871024,75871387,0]], stop_blocks=[], cds_sha256='exact-CDS',
            identification=dict(proposed_name='NANOG', near_best=[dict(gene_symbol='NANOG')]))
        result.update(changes)
        return result

    def test_named_chicken_duplicate_is_not_relabelled_by_nanog_hit(self):
        record = self.chicken(); before = deepcopy(record); found = identity(record)
        self.assertEqual(found['label'], 'NANOG-like')
        self.assertEqual(found['family_hint'], 'NANOG')
        self.assertTrue(found['database_name_differs'])
        self.assertFalse(found['orthology_confirmed'])
        self.assertEqual(found['function_status'], 'not_assessed')
        self.assertEqual(record, before)

    def test_species_scope_and_source_ids_are_required(self):
        record = self.chicken(genome_id='Homo_sapiens_GCF_000001405.40')
        self.assertEqual(identity(record)['status'], 'homology_only')
        record = self.chicken(annotation_name='LOC999', original_gene='gene-LOC999')
        self.assertEqual(identity(record)['status'], 'homology_only')

    def test_query_species_list_is_not_a_target_species(self):
        record=self.chicken(genome_id='custom_target',species=['Gallus gallus'])
        self.assertEqual(identity(record)['status'],'homology_only')
        record['species']='Gallus gallus'
        self.assertEqual(identity(record)['status'],'curated_locus_label')

    def test_named_paralogs_and_like_suffixes_are_preserved(self):
        for label in ['NANOG', 'NANOG-like', 'POU5F3', 'POU5F1', 'MYCL', 'MYC']:
            found = identity(self.chicken(annotation_name=label, original_gene='gene-'+label))
            self.assertEqual(found['label'], label)
            self.assertEqual(found['basis'], 'provided_annotation')

    def test_unknown_loc_with_a_strong_hit_is_still_a_family_candidate(self):
        found = identity(self.chicken(annotation_name='LOC999', original_gene='gene-LOC999'))
        self.assertEqual(found['display_name'], 'NANOG family candidate (LOC999)')
        self.assertEqual(found['status'], 'homology_only')

    def test_accession_without_a_symbol_does_not_invent_a_family(self):
        found = identity(self.chicken(annotation_name='LOC999', original_gene='gene-LOC999',
            identification=dict(proposed_name='P12345', near_best=[dict(accession='P12345',gene_symbol='')])))
        self.assertEqual(found['status'], 'unresolved')
        self.assertEqual(found['label'], 'LOC999')

    def test_exact_prediction_can_retain_annotated_locus_identity(self):
        annotation = self.chicken()
        prediction = self.chicken(model_id='prediction', source='tiberius', annotation_name='', original_gene='g1')
        data = dict(models=[annotation, prediction]); before = deepcopy(data)
        found = model_identities(data)
        self.assertEqual(found['prediction']['label'], 'NANOG-like')
        self.assertEqual(found['prediction']['basis'], 'exact_native_cds_and_sequence')
        self.assertEqual(data, before)

    def test_overlap_or_equal_protein_cannot_transfer_a_copy_name(self):
        annotation = self.chicken()
        for changes in [dict(cds_sha256='other-CDS'), dict(cds_sha256=''),
                        dict(strand='-'), dict(contig='another-contig'),
                        dict(cds_blocks=[[75862830,75862914,1],[75870577,75870700,0],[75871024,75871387,0]]),
                        dict(stop_blocks=[[75871387,75871390]])]:
            prediction = self.chicken(model_id='prediction', source='tiberius', annotation_name='', original_gene='g1', **changes)
            with self.subTest(changes=changes):
                found = model_identities(dict(models=[annotation,prediction]))['prediction']
                self.assertEqual(found['status'], 'homology_only')

    def test_identical_annotations_with_conflicting_gene_ids_stay_ambiguous(self):
        a = self.chicken(); b = self.chicken(model_id='other_gene', original_gene='gene-NANOG', annotation_name='NANOG')
        p = self.chicken(model_id='prediction', source='tiberius', original_gene='g1', annotation_name='')
        found = model_identities(dict(models=[a,b,p]))['prediction']
        self.assertEqual(found['status'], 'homology_only')


if __name__ == '__main__':
    unittest.main()
