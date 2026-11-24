"""Native Stage4 prepare/apply preserves sources and the exact review gate.

The old native fixture is extended with synthetic taxa. Tree construction is
stubbed at its engine boundary; real native dossier/prepare/apply code runs.
"""
from contextlib import redirect_stdout
from copy import deepcopy
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.sequence import read_fasta, write_fasta, write_tsv, sha256
from genecontext.stage3_analysis import read_table
from genecontext.stage4 import run
from genecontext.stage4_phylogeny import unique_alignment
from genecontext.tree_review_apply import prepare, apply
from genecontext.tree_review_service import default_review, validate_panel, prepare_overlay, ReviewState
from tests.test_stage4 import fixture, seal, digest


def four_taxa_source(root):
    """Extend the existing two-model Stage3 contract fixture without shortcuts."""
    state, proteins = fixture(root, reviewed=True)
    coding = read_fasta(root / 'selected.cds.fna')
    catalogue = read_table(root / 'sequence_catalogue.tsv')
    quality = read_table(root / 'model_quality.tsv'); references = read_table(root / 'reference_comparisons.tsv')
    comparison = json.loads((root / 'comparison.json').read_text())
    evidence = json.loads((root / 'finalization.evidence.json').read_text())
    for mid, amino, codon in [('mC', 'G', 'GGT'), ('mD', 'C', 'TGT')]:
        proteins[mid] = 'MA' + amino * 35; coding[mid] = 'ATGGCT' + codon * 35 + 'TAA'
        row = dict(catalogue[0], model_id=mid, genome_id=mid + '_GCF_1.1', physical_id=mid,
            species=mid[-1], original_model='rna_' + mid, protein_sha256=digest(proteins[mid]), cds_sha256=digest(coding[mid]),
            protein_file='sequences/' + mid + '.protein.faa', cds_file='sequences/' + mid + '.cds.fna')
        catalogue.append(row); quality.append(dict(quality[0], model_id=mid)); references.append(dict(references[0], model_id=mid))
        comparison['selected'].append(dict(model_id=mid, group='g', genome_id=row['genome_id']))
        evidence['models'][mid] = deepcopy(evidence['models']['mA'])
        evidence['spec']['buckets'].append(dict(key=mid, gene='TEST', group='g', genome_id=row['genome_id'], species=mid[-1], model_ids=[mid]))
        state['choices'][mid] = dict(primary=mid, reviewed=True, reason='Explicit synthetic choice')
    for row in catalogue: row['species'] = row['model_id'][-1]
    for bucket in evidence['spec']['buckets']: bucket['species'] = bucket['key'][-1]
    signature = [('g', r['genome_id'], r['model_id'], r['cds_sha256']) for r in catalogue]
    evidence['spec']['set_id'] = digest(json.dumps(['r', sorted(signature)], sort_keys=True))
    state['set_id'] = evidence['spec']['set_id']
    for name, value in [('comparison.json', comparison), ('finalization.evidence.json', evidence), ('finalization.decisions.json', state)]:
        (root / name).write_text(json.dumps(value))
    for name, rows in [('sequence_catalogue.tsv', catalogue), ('model_quality.tsv', quality), ('reference_comparisons.tsv', references)]:
        write_tsv(root / name, rows, list(rows[0]))
    write_fasta(root / 'selected.proteins.faa', proteins.items()); write_fasta(root / 'selected.cds.fna', coding.items())
    write_fasta(root / 'queries/g/proteins.aligned.faa', list(proteins.items()) + [('reference_q', proteins['mA'])])
    seal(root)


def synthetic_tree(folder, aligned, labels, *args):
    masked, leaves, mask = unique_alignment(aligned, labels)
    alias = {mid: leaf['tree_id'] for leaf in leaves for mid in leaf['sequence_ids']}
    text = '((' + alias['mA'] + ':1,' + alias['mB'] + ':1)80:1,(' + alias['mC'] + ':1,' + alias['mD'] + ':1)80:1);'
    write_fasta(folder / 'tree.aligned.faa', masked.items()); (folder / 'ml.treefile').write_text(text)
    return dict(status='computed', leaves=leaves, mask=mask, rooting='unrooted', root_basis='',
        outgroup_sequence_ids=[], bootstrap_replicates=1000, method='Synthetic fixed tree fixture', commands=[])


class NativeReviewApply(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.stage3 = self.root / 'stage3'; self.stage3.mkdir(); four_taxa_source(self.stage3)
        self.reference = self.root / 'species.nwk'; self.reference.write_text('((A,B),(C,D));')
        self.native = self.root / 'native'
        args = parser().parse_args(['stage4', '--results', str(self.stage3), '--out', str(self.native), '--species-tree', str(self.reference)])
        with patch('genecontext.stage4.build_tree', side_effect=synthetic_tree), redirect_stdout(StringIO()): run(args)
        self.native_hashes = self.sealed_hashes(self.native)
        self.stage3_hashes = self.sealed_hashes(self.stage3)

    def tearDown(self): self.tmp.cleanup()

    def sealed_hashes(self, root):
        return {n: sha256(root / n) for n in json.loads((root / 'seal.json').read_text())} | {'seal.json': sha256(root / 'seal.json')}

    def prepared(self):
        output = self.root / 'review-copy'
        prepare(SimpleNamespace(results=str(self.native), out=str(output), species_tree=str(self.reference), species_lineages=None))
        inputs = json.loads((output / 'tree_inputs.json').read_text())
        return output, inputs, default_review(inputs)

    def apply_review(self, prepared, review, name='applied'):
        file = self.root / (name + '-review.json'); file.write_text(json.dumps(review))
        output = self.root / name
        apply(SimpleNamespace(results=str(prepared), reviews=str(file), out=str(output)))
        return output, json.loads((output / 'evidence.json').read_text())

    def approve(self, inputs, review):
        for key, panel in review['panels'].items():
            preview = validate_panel(inputs, panel)
            for gid, group in panel['groups'].items():
                group.update(approved=True, confirmed_signature=preview['groups'][gid]['signature'])

    def test_native_prepare_uses_all_four_real_model_ids_and_source_is_unchanged(self):
        prepared, inputs, review = self.prepared()
        self.assertEqual(set(next(iter(inputs['panels']['analysis']['groups'].values()))['models'][i]['taxon'] for i in range(4)), set('ABCD'))
        result = prepare_overlay(review, inputs)['panels']['analysis']
        self.assertEqual({r['proposed_points'] for r in result['models'].values()}, {12.0})
        self.assertTrue(all(r['points'] is None for r in result['models'].values()))
        self.assertEqual(self.sealed_hashes(self.native), self.native_hashes)
        self.assertEqual(self.sealed_hashes(self.stage3), self.stage3_hashes)

    def test_approved_apply_enables_t_preserves_sources_and_seed_survives(self):
        prepared, inputs, review = self.prepared(); self.approve(inputs, review)
        before = self.sealed_hashes(prepared)
        output, data = self.apply_review(prepared, review)
        self.assertEqual({r['tree_score'] for r in data['rows']}, {12.0})
        self.assertTrue(all(r['tree_reviewed'] for r in data['rows']))
        state = ReviewState(output)
        restored = prepare_overlay(state.reviews(), state.manifest)['panels']['analysis']
        self.assertEqual({r['points'] for r in restored['models'].values()}, {12.0})
        self.assertEqual(self.sealed_hashes(prepared), before)
        self.assertEqual(self.sealed_hashes(self.native), self.native_hashes)
        self.assertEqual(self.sealed_hashes(self.stage3), self.stage3_hashes)
        self.assertEqual(read_fasta(output / 'primary_targets.proteins.faa'), read_fasta(prepared / 'primary_targets.proteins.faa'))

    def test_pending_apply_keeps_t_na_and_preserves_proposed(self):
        prepared, inputs, review = self.prepared()
        output, data = self.apply_review(prepared, review)
        self.assertTrue(all(r['tree_score'] is None and r['tree_proposed_score'] == 12 for r in data['rows']))
        self.assertTrue(all(not r['tree_reviewed'] for r in data['rows']))
        state = ReviewState(output)
        self.assertTrue(all(r['points'] is None for r in prepare_overlay(state.reviews(), state.manifest)['panels']['analysis']['models'].values()))
        with (output / 'final_table.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertTrue(all(r['tree_score'] == 'NA' for r in rows))

    def test_changed_outgroup_apply_clears_confirmation_and_keeps_raw_t(self):
        prepared, inputs, review = self.prepared(); self.approve(inputs, review)
        group = next(iter(review['panels']['analysis']['groups'].values()))
        leaf = next(iter(inputs['panels']['analysis']['groups'].values()))['models'][0]['tree_id']
        group.update(outgroup_ids=[leaf], root_basis='Synthetic orientation check')
        output, data = self.apply_review(prepared, review)
        self.assertTrue(all(r['tree_score'] is None and r['tree_proposed_score'] == 12 for r in data['rows']))
        state = ReviewState(output)
        self.assertFalse(next(iter(prepare_overlay(state.reviews(), state.manifest)['panels']['analysis']['groups'].values()))['approved'])

    def test_client_numbers_are_not_copied_into_native_dossier(self):
        prepared, inputs, review = self.prepared(); self.approve(inputs, review)
        for group in review['panels']['analysis']['groups'].values():
            group['points'] = 15; group['placement'] = {'models': {'mA': {'points': 15}}}
        _, data = self.apply_review(prepared, review)
        self.assertEqual({r['tree_score'] for r in data['rows']}, {12.0})

    def test_tampered_source_is_rejected_without_creating_output(self):
        prepared, inputs, review = self.prepared(); (prepared / 'report.html').write_text('tampered')
        with self.assertRaisesRegex(ValueError, 'seal differs'): self.apply_review(prepared, review)
        self.assertFalse((self.root / 'applied').exists())


if __name__ == '__main__': unittest.main()
