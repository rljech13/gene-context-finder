import unittest

from genecontext.sequence import Model
from genecontext.stage2 import model_record
from genecontext.stage2_intervals import canonical_windows, project_model


class SharedPredictionIntervals(unittest.TestCase):
    def regions(self):
        return [dict(region_id='a',locus_id='a',genome_id='g',gene='TEST',contig='chr',
                     strand='+',window_start=100,window_end=1100,anchor_start=500,anchor_end=600),
                dict(region_id='b',locus_id='b',genome_id='g',gene='TEST',contig='chr',
                     strand='-',window_start=500,window_end=1500,anchor_start=800,anchor_end=900)]

    def test_order_independent_shared_input_and_projection(self):
        regions=self.regions();intervals=canonical_windows(regions,2000)
        self.assertEqual(intervals,canonical_windows(list(reversed(regions)),2000))
        self.assertEqual(len(intervals),1);c=intervals[0]
        self.assertEqual((c['window_start'],c['window_end']),(100,1500))
        model=Model('p',c['region_id'],'-',cds=[(700,790,0),(450,510,0)],gene='p')
        records=[model_record(project_model(model,c,r),r,'tiberius','ATG','M','partial',[]) for r in regions]
        self.assertEqual(records[0]['cds_blocks'],records[1]['cds_blocks'])
        self.assertEqual(records[0]['gene_group'],records[1]['gene_group'])

    def test_transitive_overlap_is_bounded_not_silently_split(self):
        rows=self.regions()
        rows.append(dict(rows[1],region_id='c',locus_id='c',window_start=1400,window_end=2400))
        with self.assertRaisesRegex(ValueError,'Overlapping prediction windows'):
            canonical_windows(rows,2000)

    def test_separate_contigs_and_genomes_never_merge(self):
        rows=self.regions();rows.extend([dict(rows[0],region_id='c',genome_id='other'),
                                        dict(rows[0],region_id='d',contig='other')])
        self.assertEqual(len(canonical_windows(rows,2000)),3)

    def test_full_cds_is_preserved_at_display_edge_and_outside_is_omitted(self):
        rows=self.regions();c=canonical_windows(rows,2000)[0]
        model=Model('p',c['region_id'],'+',cds=[(300,420,0)],gene='p')
        projected=project_model(model,c,rows[1]);self.assertLess(projected.start,0)
        rec=model_record(projected,rows[1],'tiberius','ATG','M','partial',[])
        self.assertEqual(rec['start'],400);self.assertTrue(rec['window_edge'])
        outside=Model('q',c['region_id'],'+',cds=[(10,100,0)])
        self.assertIsNone(project_model(outside,c,rows[1]))
