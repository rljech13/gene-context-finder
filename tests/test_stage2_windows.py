"""The preliminary neighborhood defines real extraction, not a display stretch."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.sequence import write_tsv
from genecontext.stage2 import extract_intervals, plan
from genecontext.stage2_windows import resolve_windows
from genecontext.stage2_intervals import canonical_windows


class NeighborhoodBounds(unittest.TestCase):
    def fixture(self):
        genes = [dict(gene_id='left',name='LEFT',contig='chr',start=100,end=250,strand='-',side='left'),
                 dict(gene_id='right',name='RIGHT',contig='chr',start=800,end=950,strand='+',side='right')]
        locus = dict(genome_id='g',locus_id='l',contig='chr',start=400,end=600,contig_length=2000,
                     early_synteny=dict(genome_id='g',locus_id='l',window_start=100,window_end=950,
                                        neighbors=genes,display_only_genes=[],status='two_sided_context_support'))
        region = dict(genome_id='g',locus_id='l',region_id='l',contig='chr',strand='-',
                      anchor_start=410,anchor_end=590,window_start=360,window_end=640)
        return dict(loci=[locus]), [region], [dict(genome_id='g',locus_id='l')]

    def test_full_outer_genes_in_native_coordinates_and_input_unchanged(self):
        data,regions,rows=self.fixture();before=deepcopy((data,regions,rows))
        out=resolve_windows(data,regions,rows)
        self.assertEqual((out[0]['window_start'],out[0]['window_end']),(100,950))
        self.assertEqual(out[0]['window_source'],'pre-synteny')
        self.assertEqual(len(out[0]['window_resolution']['neighbor_genes']),2)
        self.assertEqual((data,regions,rows),before)
        self.assertEqual(canonical_windows(out,2000)[0]['window_end'],950)

    def test_display_only_gene_and_reviewed_anchor_are_not_cut(self):
        data,regions,rows=self.fixture()
        data['loci'][0]['early_synteny']['display_only_genes']=[dict(
            gene_id='overlap',contig='chr',start=50,end=650,strand='+')]
        regions[0].update(anchor_start=410,anchor_end=1000,window_end=1050)
        out=resolve_windows(data,regions,rows)[0]
        self.assertEqual((out['window_start'],out['window_end']),(50,1000))
        self.assertEqual(out['window_resolution']['display_only_genes'][0]['gene_id'],'overlap')

    def test_manual_bounds_and_explicit_legacy_policy_take_precedence(self):
        data,regions,rows=self.fixture()
        rows[0].update(window_start='360',window_end='640')
        out=resolve_windows(data,regions,rows)[0]
        self.assertEqual((out['window_start'],out['window_end'],out['window_source']),(360,640,'manual'))
        rows[0].clear();rows[0].update(genome_id='g',locus_id='l')
        self.assertEqual(resolve_windows(data,regions,rows,'fixed')[0]['window_source'],'fixed')

    def test_absent_annotation_uses_explicit_fallback_not_target_only_window(self):
        data,regions,rows=self.fixture()
        data['loci'][0]['early_synteny'].update(neighbors=[],status='annotation_unavailable',window_start=400,window_end=600)
        out=resolve_windows(data,regions,rows)[0]
        self.assertEqual((out['window_start'],out['window_end']),(360,640))
        self.assertEqual(out['window_source'],'fixed_fallback')
        self.assertEqual(out['window_resolution']['preliminary_status'],'annotation_unavailable')

    def test_stale_identity_wrong_contig_and_outside_gene_are_errors(self):
        for mutation in ['identity','contig','end']:
            data,regions,rows=self.fixture();c=data['loci'][0]['early_synteny']
            if mutation=='identity':c['locus_id']='other'
            elif mutation=='contig':c['neighbors'][0]['contig']='other'
            else:c['neighbors'][0]['end']=2001
            with self.subTest(mutation=mutation),self.assertRaisesRegex(ValueError,'Preliminary|preliminary'):
                resolve_windows(data,regions,rows)

    def test_resolved_neighborhood_is_subject_to_prediction_cap(self):
        data,regions,rows=self.fixture()
        with self.assertRaisesRegex(ValueError,'exceeds'):
            canonical_windows(resolve_windows(data,regions,rows),500)

    def test_cli_plan_needs_no_genome_or_predictor_and_bounds_drive_extraction(self):
        data,_,_=self.fixture();data.update(review_id='review',config=dict(gene='X',window=50),
                                           genomes=[dict(genome_id='g',status='success')])
        data['loci'][0].update(gene='X',strand='-',pass_filters=True,boundary_source='annotation')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'stage1';source.mkdir()
            (source/'evidence.json').write_text(json.dumps(data))
            write_tsv(root/'decisions.tsv',[dict(review_id='review',genome_id='g',locus_id='l',decision='inspect')],
                      ['review_id','genome_id','locus_id','decision'])
            args=parser().parse_args(['stage2-plan','--results',str(source),'--decisions',str(root/'decisions.tsv'),
                                       '--out',str(root/'plan')])
            with patch('genecontext.stage2.verify',return_value=(source,data)):
                self.assertEqual(plan(args),0)
            saved=json.loads((root/'plan/window_plan.json').read_text());r=saved['regions'][0]
            self.assertEqual((r['window_start'],r['window_end']),(100,950))
            self.assertEqual(saved['total_prediction_bases'],850)
            self.assertEqual({p.name for p in (root/'plan').iterdir()},{'window_plan.json','windows.tsv'})
            (root/'genome.fa').write_text('>chr\n'+'A'*100+'C'*850+'G'*1050+'\n')
            dna,_=extract_intervals(root/'genome.fa',{r['region_id']:(r['contig'],r['window_start'],r['window_end'])})
            self.assertEqual(dna['l'],'C'*850)


if __name__=='__main__':unittest.main()
