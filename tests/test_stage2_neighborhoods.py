"""Fresh coding-gene coverage, and refusal before any DNA/predictor work."""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.sequence import sha256, write_tsv
from genecontext.stage2 import plan, planned_regions, run
from genecontext.stage2_neighborhoods import annotation_neighborhoods
from genecontext.stage2_windows import resolve_windows


def fixture(root, count=10):
    gff=root/'annotation.gff3'
    records=[('target',5_000_000,5_010_000)]+[(f'L{i}',5_000_000-i*300_000,5_003_000-i*300_000) for i in range(1,count+1)]+[(f'R{i}',5_010_000+i*300_000,5_013_000+i*300_000) for i in range(1,count+1)]
    lines=['##gff-version 3']
    for name,a,b in records:
        lines.extend([f'chr\tx\tgene\t{a+1}\t{b}\t.\t+\t.\tID={name};Name={name};gene_biotype=protein_coding',
                      f'chr\tx\tmRNA\t{a+1}\t{b}\t.\t+\t.\tID={name}.t;Parent={name}',
                      f'chr\tx\tCDS\t{a+1}\t{a+301}\t.\t+\t0\tParent={name}.t'])
    # An extra isoform and a nearby noncoding gene must not fill neighbor slots.
    lines.extend(['chr\tx\tmRNA\t4700001\t4703000\t.\t+\t.\tID=L1.other;Parent=L1',
                  'chr\tx\tCDS\t4700101\t4700300\t.\t+\t0\tParent=L1.other',
                  'chr\tx\tgene\t4990001\t4991000\t.\t+\t.\tID=RNA;gene_biotype=lncRNA'])
    gff.write_text('\n'.join(lines)+'\n')
    locus=dict(genome_id='g',locus_id='l',gene='FAMILY',contig='chr',start=5_000_000,end=5_010_000,strand='+',contig_length=10_000_000,annotation_gene_id='target',pass_filters=True,boundary_source='annotation',
               early_synteny=dict(genome_id='g',locus_id='l',status='insufficient_anchors',window_start=4_700_000,window_end=5_010_000,
                neighbors=[dict(gene_id='L1',name='L1',contig='chr',start=4_700_000,end=4_703_000,strand='+',side='left')],display_only_genes=[]))
    data=dict(review_id='review',config=dict(gene='FAMILY',window=50_000),genomes=[dict(genome_id='g',status='completed',annotation_inputs=dict(gff3=dict(sha256=sha256(gff))))],loci=[locus])
    region=dict(genome_id='g',locus_id='l',region_id='l',contig='chr',strand='+',anchor_start=locus['start'],anchor_end=locus['end'],window_start=4_950_000,window_end=5_060_000)
    rows=[dict(genome_id='g',locus_id='l')]
    manifest=root/'genomes.tsv';write_tsv(manifest,[dict(genome_id='g',fasta='not-read.fna',gff3=gff.name)],['genome_id','fasta','gff3'])
    source=root/'source';source.mkdir();(source/'evidence.json').write_text(json.dumps(data))
    decisions=root/'choices.tsv';write_tsv(decisions,[dict(review_id='review',genome_id='g',locus_id='l',decision='inspect',reason='check_context')],['review_id','genome_id','locus_id','decision','reason'])
    return data,[region],rows,manifest,source,decisions


class FullNeighborhoods(unittest.TestCase):
    def test_sparse_old_map_expands_to_ten_distinct_coding_genes_per_side(self):
        with TemporaryDirectory() as tmp:
            data,regions,rows,manifest,_,_=fixture(Path(tmp));before=deepcopy((data,regions,rows))
            contexts=annotation_neighborhoods(data,regions,manifest)
            out=resolve_windows(data,regions,rows,'annotation-neighbors',contexts)[0]
            self.assertEqual((out['window_start'],out['window_end']),(2_000_000,8_013_000))
            coverage=out['window_resolution']['neighborhood_coverage']
            self.assertEqual(coverage['captured'],dict(left=10,right=10))
            self.assertTrue(coverage['ready_for_prediction'])
            self.assertEqual(out['window_resolution']['previous_preliminary_neighbors'],1)
            self.assertNotIn('RNA',coverage['fully_captured_neighbor_ids'])
            self.assertEqual((data,regions,rows),before)

    def test_large_plan_retains_complete_genes_and_blocks_real_execution(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);data,_,_,manifest,source,decisions=fixture(root)
            args=parser().parse_args(['stage2-plan','--results',str(source),'--decisions',str(decisions),'--genomes',str(manifest),'--out',str(root/'plan')])
            with patch('genecontext.stage2.verify',return_value=(source,data)):self.assertEqual(plan(args),0)
            saved=json.loads((root/'plan/window_plan.json').read_text())
            self.assertFalse(saved['ready_for_prediction']);self.assertEqual(saved['regions'][0]['window_end'],8_013_000)
            self.assertEqual(saved['execution_blockers'][0]['kind'],'window_size')
            with self.assertRaisesRegex(ValueError,'exceeds --max-window-bases'):planned_regions(source,data,args)
            self.assertFalse((root/'not-read.fna').exists())

    def test_missing_annotation_is_explicit_and_predictor_is_never_called(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);data,_,_,manifest,source,decisions=fixture(root)
            write_tsv(manifest,[dict(genome_id='g',fasta='not-read.fna',gff3='')],['genome_id','fasta','gff3'])
            args=parser().parse_args(['stage2','--results',str(source),'--decisions',str(decisions),'--genomes',str(manifest),'--protein-db','not-read-db','--out',str(root/'result')])
            with patch('genecontext.stage2.verify',return_value=(source,data)),patch('genecontext.stage2.verify_db') as db,patch('genecontext.tiberius_backend.predict_batches') as predictor:
                with self.assertRaisesRegex(ValueError,'annotation_unavailable'):run(args)
                db.assert_not_called();predictor.assert_not_called()
            self.assertFalse((root/'result').exists())

    def test_manual_window_is_preserved_but_missing_neighbors_block_prediction(self):
        with TemporaryDirectory() as tmp:
            data,regions,rows,manifest,_,_=fixture(Path(tmp));rows[0].update(window_start='4950000',window_end='5060000')
            contexts=annotation_neighborhoods(data,regions,manifest)
            out=resolve_windows(data,regions,rows,'annotation-neighbors',contexts)[0]
            self.assertEqual((out['window_start'],out['window_end']),(4_950_000,5_060_000))
            self.assertFalse(out['window_resolution']['neighborhood_coverage']['ready_for_prediction'])
            self.assertEqual(len(out['window_resolution']['neighborhood_coverage']['neighbors_outside_reviewed_window']),20)

    def test_source_annotation_hash_is_checked_before_using_new_neighbors(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);data,regions,_,manifest,_,_=fixture(root)
            with (root/'annotation.gff3').open('a') as out:out.write('# different release\n')
            with self.assertRaisesRegex(ValueError,'differs from stage 1'):annotation_neighborhoods(data,regions,manifest)

    def test_contig_shortage_does_not_invent_neighbors(self):
        with TemporaryDirectory() as tmp:
            data,regions,rows,manifest,_,_=fixture(Path(tmp),count=3)
            contexts=annotation_neighborhoods(data,regions,manifest)
            out=resolve_windows(data,regions,rows,'annotation-neighbors',contexts)[0]
            c=out['window_resolution']['neighborhood_coverage']
            self.assertEqual(c['captured'],dict(left=3,right=3));self.assertFalse(c['ready_for_prediction'])
            self.assertEqual(c['limits'],dict(left='contig_or_annotation_end',right='contig_or_annotation_end'))

    def test_default_requires_inventory_legacy_policy_stays_explicit(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);data,_,_,_,source,decisions=fixture(root)
            args=parser().parse_args(['stage2-plan','--results',str(source),'--decisions',str(decisions),'--out',str(root/'plan')])
            with self.assertRaisesRegex(ValueError,'requires --genomes'):planned_regions(source,data,args)
            args.window_policy='pre-synteny';regions,_,_,_=planned_regions(source,data,args)
            self.assertEqual(regions[0]['window_start'],4_700_000)
