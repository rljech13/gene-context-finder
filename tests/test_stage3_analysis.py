"""Scientific boundary and interpretation checks for the final report."""
import unittest

from genecontext.stage3_analysis import extract_flanks,native_interval,pair_relation


class ModelIntervals(unittest.TestCase):
    def test_flanks_follow_selected_model_on_both_strands(self):
        row=dict(window_start=100,window_end=116,start=104,end=108,strand='+')
        seq='AAAACCCCGGGGTTTT'
        flanks,info=extract_flanks(row,seq,requested=3)
        self.assertEqual(flanks,dict(upstream='AAA',downstream='GGG'))
        self.assertEqual((info['upstream_native_start'],info['upstream_native_end']),(101,104))
        self.assertEqual((info['downstream_native_start'],info['downstream_native_end']),(108,111))
        row['strand']='-';flanks,info=extract_flanks(row,seq,requested=3)
        self.assertEqual(flanks,dict(upstream='CCC',downstream='TTT'))
        self.assertEqual((info['upstream_native_start'],info['upstream_native_end']),(108,111))
        self.assertEqual((info['downstream_native_start'],info['downstream_native_end']),(101,104))

    def test_window_edges_are_reported_without_inventing_flank_sequence(self):
        row=dict(window_start=100,window_end=116,start=104,end=108,strand='-')
        flanks,info=extract_flanks(row,'AAAACCCCGGGGTTTT',requested=50000)
        self.assertEqual((len(flanks['upstream']),len(flanks['downstream'])),(8,4))
        self.assertEqual(info['flank_status'],'limited_by_saved_window')
        row['start']=99;flanks,info=extract_flanks(row,'AAAACCCCGGGGTTTT')
        self.assertEqual(flanks,{});self.assertEqual(info['flank_status'],'model_not_fully_observed')

    def test_native_coordinates_roundtrip_in_target_orientation(self):
        row=dict(window_start=100,window_end=116,strand='+')
        self.assertEqual(native_interval(row,2,7),(102,107))
        row['strand']='-';self.assertEqual(native_interval(row,2,7),(109,114))

    def test_shared_assembly_windows_are_not_independent_conservation_evidence(self):
        a=dict(genome_id='a',contig='chr1',window_start=100,window_end=200)
        b=dict(a,window_start=180,window_end=280)
        self.assertEqual(pair_relation(a,b)['native_window_overlap_bp'],20)
        b['contig']='chr2';self.assertEqual(pair_relation(a,b)['native_window_overlap_bp'],0)
        b['contig']='chr1';b['genome_id']='b'
        self.assertEqual(pair_relation(a,b)['pair_scope'],'between_genomes')
        self.assertEqual(pair_relation(a,b)['native_window_overlap_bp'],0)
