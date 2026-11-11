import copy
from pathlib import Path
import tempfile
import unittest
from genecontext.stage2_display import automatic_comments, filtered_fragments

class DisplayEvidence(unittest.TestCase):
    def fixture(self):
        a=dict(model_id='a',region_id='r',gene_group='g',genome_id='bird',source='annotation',strand='+',overlaps_anchor=True,annotation_name='TARGET',protein_length=462,cds_blocks=[[100,805,0],[100000,100112,0],[160000,160131,2],[204000,204165,0],[214000,214276,0]])
        t=dict(a,model_id='t',source='tiberius',protein_length=191,cds_blocks=[[160000,160131,2],[204000,204165,0],[214000,214276,0]])
        return dict(models=[a,t],regions=[dict(region_id='r',genome_id='bird',contig='chr')],prediction_intervals=[dict(region_id='p',region_ids=['r'],window_start=1000)])
    def test_size_difference_is_only_a_comment_and_uses_separate_units(self):
        data=self.fixture();before=copy.deepcopy(data);notes=automatic_comments(data)
        self.assertIn('Genomic span differs',notes['r'][0]);self.assertIn('Protein length differs',notes['r'][0]);self.assertEqual(data,before)
    def test_small_changes_do_not_generate_size_comments(self):
        d=self.fixture();d['models'][1].update(protein_length=450,cds_blocks=d['models'][0]['cds_blocks']);self.assertEqual(automatic_comments(d),{})
    def test_short_gene_absolute_floor_avoids_spurious_span_warning(self):
        d=self.fixture();d['models'][0].update(cds_blocks=[[0,3000,0]],protein_length=100);d['models'][1].update(cds_blocks=[[0,500,0]],protein_length=100);self.assertEqual(automatic_comments(d),{})
    def test_raw_removed_fragment_is_display_only_and_mapped_in_native_coordinates(self):
        d=self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'context_prediction/tiberius_001';p.mkdir(parents=True)
            # Native E2 is [100000,100112); predictor window starts at 1000.
            text='p\tT\tmRNA\t99001\t99112\t.\t+\t.\tID=raw;Parent=g\np\tT\tCDS\t99001\t99112\t.\t+\t0\tParent=raw\n'
            (p/'predictions.gff3.raw.gff3').write_text(text);(p/'predictions.gff3').write_text('##gff-version 3\n')
            before=copy.deepcopy(d);out=filtered_fragments(tmp,d);self.assertEqual(len(out),1);m=out[0];self.assertEqual(m['cds_blocks'],[[100000,100112,0]]);self.assertFalse(m['reviewable']);self.assertTrue(m['diagnostic_fragment']);self.assertEqual(m['target_cds_overlap'],112);self.assertEqual(d,before)
            (p/'predictions.gff3').write_text(text);self.assertEqual(filtered_fragments(tmp,d),[])
    def test_raw_unrelated_or_missing_outputs_are_not_invented(self):
        with tempfile.TemporaryDirectory() as tmp:self.assertEqual(filtered_fragments(tmp,self.fixture()),[])

if __name__=='__main__':unittest.main()
