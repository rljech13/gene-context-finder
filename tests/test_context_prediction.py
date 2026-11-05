import argparse
import csv
import json
import random
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from genecontext.cli import parser
from genecontext.context_predict import identify_proteins, native_interval, prediction_model
from genecontext.sequence import read_fasta, read_gff, reverse_complement, translate, write_fasta
from genecontext.tools import Runner


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class PredictionUnits(unittest.TestCase):
    def test_model_override_and_native_mapping(self):
        args = argparse.Namespace(augustus_species='human')
        self.assertEqual(prediction_model(args, dict(augustus_species='chicken')), 'chicken')
        self.assertEqual(prediction_model(args, {}), 'human')
        self.assertEqual(native_interval(dict(strand='-', window_start=100, window_end=1000), 20, 70), (930, 980))
        self.assertEqual(native_interval(dict(strand='+', window_start=100, window_end=1000), 20, 70), (120, 170))


@unittest.skipUnless(all(shutil.which(p) for p in ['augustus', 'samtools', 'miniprot', 'blastn', 'blastp', 'blastx', 'makeblastdb', 'mafft']), 'requires actual prediction and comparison tools')
class PredictionIntegration(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        source = Path(__file__).parent / 'fixtures/HS04636.fna'
        self.sequence = next(iter(read_fasta(source).values()))

    def tearDown(self):
        self.temp.cleanup()

    def test_prediction_precedes_identification_with_and_without_gff(self):
        genome = self.root / 'genome.fna'; write_fasta(genome, [('chr', self.sequence)])
        reverse = self.root / 'reverse.fna'; write_fasta(reverse, [('chr', reverse_complement(self.sequence))])
        # Deliberately fictitious annotation: it must not create a predicted gene.
        gff = self.root / 'annotation.gff3'
        gff.write_text('##gff-version 3\nchr\tfixture\tgene\t8000\t8500\t.\t+\t.\tID=FICTITIOUS;Name=FICTITIOUS\n')
        inputs = self.root / 'genomes.tsv'
        inputs.write_text('genome_id\tfasta\tgff3\taugustus_species\nannotated\t'+str(genome)+'\t'+str(gff)+'\thuman\nraw\t'+str(reverse)+'\t\thuman\n')
        query = self.root / 'query.fna'; write_fasta(query, [('noncoding_seed', self.sequence[:400])])
        output = self.root / 'analysis'
        command = [sys.executable, '-m', 'genecontext', 'analyze', '--query', str(query), '--query-type', 'dna',
                   '--genomes', str(inputs), '--out', str(output), '--threads', '1', '--window', '50000', '--context-reference', 'none']
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
        root = output / 'comparison'
        self.assertEqual(json.loads((root / 'manifest.json').read_text())['status'], 'completed')
        features = table(root / 'neighbors.tsv')
        self.assertEqual({r['genome_id'] for r in features}, {'annotated', 'raw'})
        self.assertTrue(all(r['evidence'] == 'ab_initio_prediction' for r in features))
        self.assertTrue(all(r['name'] == 'uncharacterized' for r in features))
        self.assertTrue(all(r['identification_status'] == 'no_passing_match' for r in features))
        # The same oriented DNA must produce the same predictions despite GFF3 availability.
        a = [(r['relative_start'], r['relative_end'], r['relative_strand']) for r in features if r['genome_id']=='annotated']
        b = [(r['relative_start'], r['relative_end'], r['relative_strand']) for r in features if r['genome_id']=='raw']
        self.assertEqual(a, b)
        self.assertGreater(len(a), 0)
        self.assertNotIn('FICTITIOUS', {r['name'] for r in features})
        native = {'annotated': self.sequence, 'raw': reverse_complement(self.sequence)}
        exported = read_fasta(root / 'context_prediction_cds.fna')
        proteins = read_fasta(root / 'context_prediction_proteins.faa')
        for gid in native:
            models, _ = read_gff(root / 'context_predictions' / (gid+'.gff3'))
            self.assertTrue(any(len(m.cds)>1 for m in models))
            for model in models:
                seq = ''.join(native[gid][a:b] if model.strand=='+' else reverse_complement(native[gid][a:b]) for a,b,_ in model.cds)
                pid = model.id.removesuffix('.t1')
                self.assertEqual(seq, exported[pid])
                self.assertEqual(translate(seq).removesuffix('*'), proteins[pid])
        # Unnamed predictions still support protein homology / local gene context.
        self.assertTrue(all(r['homology_group'] for r in features))
        self.assertGreater(int(table(root / 'synteny.tsv')[0]['unique_anchors']), 0)
        self.assertEqual(len(table(root / 'provided_neighbors.tsv')), 1)

    def test_identification_retains_no_hits_and_distinguishes_ambiguous_names(self):
        # Independent sequence controls exercise the post-prediction lookup stage.
        sequence = 'M' + ''.join(random.Random(719).choices('ACDEFGHIKLMNPQRSTVWY', k=180))
        refs = self.root / 'refs.faa'
        write_fasta(refs, [('sp|ONE|FIRST GN=ALPHA', sequence), ('sp|TWO|SECOND GN=BETA', sequence)])
        args = parser().parse_args(['compare', '--results', 'unused', '--out', str(self.root), '--context-proteins', str(refs), '--context-reference', 'none', '--threads', '1'])
        metadata = {}
        assignments, _ = identify_proteins(args, self.root, {'match': sequence, 'unknown': 'M'+'W'*180}, {}, metadata, Runner(self.root))
        self.assertNotIn('unknown', assignments)
        self.assertEqual(assignments['match']['identification_status'], 'ambiguous_best_matches')
        self.assertEqual(assignments['match']['name'], 'uncharacterized')


if __name__ == '__main__':
    unittest.main()
