"""Offline report transport and preliminary-context safety regressions."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from genecontext.early_synteny import GeneIndex, compare_context, run_early_synteny
from genecontext.stage1_report import pack_report


def unpack(payload):
    def expand(value):
        if isinstance(value, list):
            return [expand(v) for v in value]
        if isinstance(value, dict):
            if set(value) == {'$gcf_ref'}:
                return expand(payload['objects'][value['$gcf_ref']])
            return {k: expand(v) for k, v in value.items()}
        return value
    return expand(payload['data'])


class ReportTransport(unittest.TestCase):
    def test_interning_preserves_all_isoforms_hsps_and_colliding_ids(self):
        model = dict(model_id='same_id', cds=[[100, 200, 0]])
        hsp = dict(hsp_id='same_hsp_id', qstart=0, qend=10)
        chains = [dict(annotation_model=model, hsps=[hsp]),
                  dict(annotation_model=copy.deepcopy(model), hsps=[copy.deepcopy(hsp)]),
                  dict(annotation_model=dict(model, cds=[[300, 400, 0]]),
                       hsps=[dict(hsp, qend=20)])]
        data = dict(review_id='test', loci=[dict(genome_id='a', chains=chains),
                                          dict(genome_id='b', chains=copy.deepcopy(chains))])
        before = copy.deepcopy(data)
        packed = pack_report(data)
        self.assertEqual(unpack(packed), data)
        self.assertEqual(data, before)
        self.assertLess(len(json.dumps(packed)), len(json.dumps(data)))
        self.assertEqual(len(unpack(packed)['loci'][0]['chains']), 3)

    def test_coordinates_variants_and_display_only_genes_round_trip(self):
        gene = dict(gene_id='copy', start=10, end=20, display_only=True,
                    models=['isoform1', 'isoform2'])
        data = dict(review_id='test', queries=[dict(id='q', members=[dict(query='a'), dict(query='b')])],
                    loci=[dict(early_synteny=dict(neighbors=[gene], display_only_genes=[gene]))],
                    report_assets=dict(dna_loci=['locus1']))
        self.assertEqual(unpack(pack_report(data)), data)


def gene(gid, start, end, strand='+', contig='chr'):
    return dict(gene_id=gid, name=gid, contig=contig, start=start, end=end,
                strand=strand, coding=True, protein_id='p_'+gid, models=[])


class ContextSafety(unittest.TestCase):
    def test_self_control_cannot_stop_cross_species_expansion(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            rows = [dict(genome_id=gid, fasta='', gff3='supplied', protein='supplied',
                         cds='', annotation_provenance='') for gid in ['a', 'b']]
            catalogs = {}
            proteins = {}
            for number, gid in enumerate(['a', 'b']):
                genes = [gene(gid+'_target', 100, 200), gene(gid+'_left', 50, 90),
                         gene(gid+'_right', 210, 250)]
                catalogs[gid] = dict(genome_id=gid, genes=genes, key=gid, cache=str(root),
                    cache_reused=True, elapsed_seconds=0, annotation_genes=3,
                    coding_genes=3, protein_genes=3, inputs={}, contigs={'chr': ('chr', 1000)})
                for i, g in enumerate(genes):
                    proteins[g['protein_id']] = 'M'+'ACDEFG'[number*3+i]*30
            origins = [dict(genome_id=gid, gene_id=gid+'_target', reference_id='r_'+gid,
                            queries=[gid], resolutions=[]) for gid in ['a', 'b']]

            def matching(panel, target, *args):
                mapping = {}
                for pid, seq in proteins.items():
                    if pid.startswith('p_'+target['genome_id']+'_'):
                        qid = 'P'+hashlib.sha256(seq.encode()).hexdigest()[:24]
                        mapping[qid] = dict(status='unique_best_match', best_matches=[dict(
                            target_protein=pid, bitscore=100, evalue=1e-20,
                            identity=1, coverage_query=1, coverage_subject=1)])
                return dict(matches=mapping, cache=str(root), cache_reused=True, panel_proteins=len(panel))

            data = dict(config=dict(query_type='protein'), queries=[],
                        genomes=[dict(genome_id='a')], loci=[dict(genome_id='a', locus_id='l',
                        annotation_gene_id='a_target', contig='chr', start=100, end=200,
                        strand='+', contig_length=1000)])
            with patch('genecontext.early_synteny.datasets', return_value=rows), \
                 patch('genecontext.early_synteny.build_catalog', side_effect=lambda row,*a: catalogs[row['genome_id']]), \
                 patch('genecontext.early_synteny.resolve_origins', return_value=(origins, [])), \
                 patch('genecontext.early_synteny.fasta_records', side_effect=lambda path: iter(proteins.items())), \
                 patch('genecontext.early_synteny.tool_versions', return_value={k:dict(path='mock') for k in ['blastp','makeblastdb']}), \
                 patch('genecontext.early_synteny.match_panel', side_effect=matching):
                result = run_early_synteny(data, genomes_manifest='unused',
                    cache_dir=root/'cache', out_dir=root/'result')
            locus = result['loci'][0]
            self.assertEqual(locus['neighbors_used'], 10)
            self.assertEqual(locus['status'], 'insufficient_anchors')
            self.assertEqual(locus['comparisons'][0]['status'], 'two_sided_context_support')
            self.assertEqual(locus['comparisons'][1]['shared_anchors'], 0)
            from genecontext.comparison import table
            self.assertEqual(table(root/'result/summary.tsv')[0]['best_shared_anchors'], '0')

    def test_other_target_and_overlapping_annotations_are_visible_not_anchors(self):
        catalog = dict(genes=[gene('anchor', 100, 200), gene('copy', 210, 220),
                              gene('nested', 130, 150, '-'), gene('left', 60, 90),
                              gene('right', 230, 250), gene('far', 9000, 9100)])
        ctx = GeneIndex(catalog).neighbors('chr', 100, 200, '+', 1, 1000,
                                           {'anchor', 'copy'}, 'anchor')
        self.assertEqual({g['gene_id'] for g in ctx['neighbors']}, {'left', 'right'})
        extras = {g['gene_id']: g for g in ctx['display_only_genes']}
        self.assertEqual(set(extras), {'copy', 'nested'})
        self.assertEqual(extras['copy']['display_reason'], 'other_query_source')
        self.assertEqual(extras['nested']['display_reason'], 'overlapping_annotation')
        self.assertTrue(all(g['display_only'] for g in extras.values()))
        reference = dict(reference_id='r', strand='+', _target_strand='+')
        c = compare_context(reference, ctx['neighbors'], ctx['neighbors'], {})
        self.assertEqual(c['shared_anchors'], 0)

    def test_reverse_orientation_preserves_neighbor_sides_and_budget(self):
        genes = [gene('a', 1000, 1100, '-'), gene('copy', 1200, 1250),
                 *[gene('left'+str(i), 100+i*50, 120+i*50) for i in range(12)],
                 *[gene('right'+str(i), 1300+i*50, 1320+i*50) for i in range(12)]]
        ctx = GeneIndex(dict(genes=genes)).neighbors('chr', 1000, 1100, '-', 10, 10000,
                                                     {'a', 'copy'}, 'a')
        for side in ['upstream', 'downstream']:
            self.assertEqual(sum(g['target_side']==side for g in ctx['neighbors']), 10)
        self.assertEqual(ctx['display_only_genes'][0]['target_side'], 'upstream')

    def test_hard_limit_rejects_eleven_before_catalog_or_blast(self):
        with tempfile.TemporaryDirectory() as d:
            with patch('genecontext.early_synteny.datasets', side_effect=AssertionError('No IO')):
                with self.assertRaisesRegex(ValueError, '<= 10'):
                    run_early_synteny({}, genomes_manifest='unused', cache_dir=Path(d)/'cache',
                                       out_dir=Path(d)/'result', config=dict(max_neighbors=11))
