"""Offline reference-panel checks: coverage, competitors and real BLASTP."""
import gzip
import json
from pathlib import Path
import random
import shutil
import tempfile
import unittest

from genecontext.cli import main, parser
from genecontext.comparison import table
from genecontext.sequence import sha256, write_fasta, write_tsv
from genecontext.target_check import check_targets, labels, proteins, summarize


class TargetCheckUnits(unittest.TestCase):
    def test_coverage_does_not_double_count_overlapping_hsps(self):
        base = dict(identity=.8, length=40, qstart=0, qend=40, sstart=100, send=140, evalue=1e-8, bitscore=80)
        repeated = dict(base, qstart=20, qend=60, sstart=120, send=160, bitscore=70)
        subject_repeated = dict(base, length=20, qstart=40, qend=60, sstart=110, send=130, bitscore=60)
        separate = dict(base, identity=.6, length=20, qstart=60, qend=80, sstart=160, send=180, bitscore=30)
        summary = summarize([repeated, subject_repeated, separate, base], 100, 200)
        self.assertEqual(summary['hsps_used'], 2)
        self.assertAlmostEqual(summary['candidate_coverage'], .6)
        self.assertAlmostEqual(summary['reference_coverage'], .3)
        self.assertAlmostEqual(summary['identity'], (40*.8+20*.6)/60)
        self.assertEqual(summary['bitscore'], 110)

    def test_labels_require_unique_exact_complete_reference_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'labels.tsv';metadata = {'R1': dict(id='a'), 'R2': dict(id='b')}
            for rows in [[dict(reference='a', family='A')],
                         [dict(reference='a', family='A'), dict(reference='a', family='A')],
                         [dict(reference='a', family='A'), dict(reference='other', family='B')]]:
                write_tsv(p, rows, ['reference', 'family'])
                with self.assertRaises(ValueError):labels(p, metadata)
            write_tsv(p, [dict(reference='a', family='A'), dict(reference='b', family='A')], ['reference', 'family'])
            self.assertEqual(labels(p, metadata), {'a': 'A', 'b': 'A'})

    def test_input_ids_stops_gzip_and_duplicate_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'p.faa.gz'
            with gzip.open(p, 'wt') as h:h.write('>a|one GN=TEST\nMACD*EFG*\n')
            seq, metadata = proteins(p, 'Q')
            self.assertEqual(seq, [('Q0000001', 'MACDXEFG')])
            self.assertEqual(metadata['Q0000001']['internal_stops'], 1)
            with self.assertRaisesRegex(ValueError, 'internal stops'):proteins(p, 'R', reference=True)
            p = Path(tmp)/'duplicate.faa';p.write_text('>a\nMACDEFG\n>a\nMACDEFG\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate'):proteins(p, 'Q')


@unittest.skipUnless(shutil.which('blastp') and shutil.which('makeblastdb'), 'BLAST+ is required')
class TargetCheckIntegration(unittest.TestCase):
    def test_real_blast_competitors_domain_matches_no_hits_and_immutable_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp);rng = random.Random(62026);alphabet = 'ACDEFGHIKLMNPQRSTVWY'
            protein = ''.join(rng.choice(alphabet) for _ in range(300))
            competitor = protein[:150]+('A' if protein[150]!='A' else 'G')+protein[151:]
            query, refs, names = root/'queries.faa', root/'refs.faa', root/'families.tsv'
            write_fasta(query, [('full|model', protein), ('domain_only', protein[110:190]),
                                ('negative', 'W'*90), ('disrupted', protein[:150]+'*'+protein[151:])])
            write_fasta(refs, [('target GN=TARGET', protein), ('same_family GN=TARGET', protein),
                               ('competitor GN=PARALOG', competitor)])
            write_tsv(names, [dict(reference='target', family='target'), dict(reference='same_family', family='target'),
                              dict(reference='competitor', family='paralog')], ['reference', 'family'])
            before = {p: sha256(p) for p in [query, refs, names]}
            args = parser().parse_args(['check-targets', '--proteins', str(query), '--references', str(refs),
                '--reference-labels', str(names), '--out', str(root/'audit'), '--threads', '1'])
            self.assertEqual(check_targets(args), 0)
            rows = {r['candidate_id']: r for r in table(root/'audit/summary.tsv')}
            self.assertEqual(set(rows), {'full|model', 'domain_only', 'negative', 'disrupted'})
            self.assertEqual(rows['full|model']['best_family'], 'target')
            self.assertEqual(rows['full|model']['alternative_family'], 'paralog')
            self.assertIn('close_other_family', rows['full|model']['flags'])
            self.assertIn('limited_coverage', rows['domain_only']['flags'])
            self.assertGreater(float(rows['domain_only']['best_candidate_coverage']), .9)
            self.assertLess(float(rows['domain_only']['best_reference_coverage']), .3)
            self.assertEqual(rows['negative']['flags'], 'no_reference_hit')
            self.assertIn('internal_stop', rows['disrupted']['flags'])
            receipt = json.loads((root/'audit/manifest.json').read_text())
            self.assertEqual(receipt['status'], 'completed')
            self.assertEqual([c['argv'][0] for c in receipt['commands']], ['makeblastdb', 'blastp'])
            self.assertIn('-db', receipt['commands'][1]['argv'])
            self.assertNotIn('-subject', receipt['commands'][1]['argv'])
            for name, digest in receipt['outputs_sha256'].items():self.assertEqual(sha256(root/'audit'/name), digest)
            for path, digest in before.items():self.assertEqual(sha256(path), digest)
            # No output overwrite and no input clobbering.
            self.assertEqual(main(['check-targets', '--proteins', str(query), '--references', str(refs),
                                   '--out', str(root/'audit')]), 2)
            self.assertEqual(main(['check-targets', '--proteins', str(query), '--references', str(refs),
                                   '--out', str(root)]), 2)

    def test_unlabelled_panel_has_no_inferred_family_or_family_confidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp);rng = random.Random(52)
            p = ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(150))
            write_fasta(root/'q.faa', [('q', p)])
            write_fasta(root/'r.faa', [('r GN=KNOWN', p)])
            self.assertEqual(main(['check-targets', '--proteins', str(root/'q.faa'), '--references', str(root/'r.faa'),
                                   '--out', str(root/'audit')]), 0)
            row = table(root/'audit/summary.tsv')[0]
            self.assertEqual(row['best_reference_name'], 'KNOWN')
            self.assertEqual(row['best_family'], '')
            self.assertEqual(row['alternative_score_ratio'], '')


if __name__ == '__main__':
    unittest.main()
