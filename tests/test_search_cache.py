"""Shared-index reuse, failure isolation and actual cold/warm search equivalence."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
from pathlib import Path
import random
import shutil
import tempfile
import unittest

from genecontext.search_cache import prepare_blast_db, prepare_genome
from genecontext.tools import Runner


VERSIONS = {'makeblastdb': {'path': '/fixture/makeblastdb', 'version': 'fixture-1'}}


class FixtureRunner:
    def __init__(self, fail=False):
        self.commands = []
        self.fail = fail

    def run(self, argv, output, log=None):
        self.commands.append({'argv': [str(x) for x in argv]})
        if self.fail:
            raise RuntimeError('Interrupted database build')
        prefix = Path(argv[argv.index('-out')+1])
        digest = hashlib.sha256(Path(argv[argv.index('-in')+1]).read_bytes()).digest()
        for suffix in ['pin', 'phr', 'psq']:
            prefix.with_suffix('.'+suffix).write_bytes(digest)
        Path(output).write_text('fixture database built\n')


class DatabaseCache(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root/'proteins.faa'
        self.fasta.write_text('>p\nMAVKLE\n')
        self.cache = self.root/'cache'

    def prepare(self, runner, **kwargs):
        return prepare_blast_db(self.fasta, self.cache, runner, 'prot',
                                versions=kwargs.pop('versions', VERSIONS), **kwargs)

    def test_same_content_at_another_path_reuses_one_database(self):
        cold = FixtureRunner()
        first, receipt = self.prepare(cold)
        self.assertFalse(receipt['reused'])
        self.assertEqual(len(cold.commands), 1)
        other = self.root/'another-query-input.faa'
        shutil.copyfile(self.fasta, other)
        warm = FixtureRunner()
        second, receipt = prepare_blast_db(other, self.cache, warm, 'prot', versions=VERSIONS)
        self.assertEqual(first, second)
        self.assertTrue(receipt['reused'])
        self.assertEqual(warm.commands, [])

    def test_input_or_tool_version_change_creates_an_independent_entry(self):
        first, _ = self.prepare(FixtureRunner())
        self.fasta.write_text('>p\nMAVKLEAA\n')
        changed, _ = self.prepare(FixtureRunner())
        upgraded, _ = self.prepare(FixtureRunner(), versions={
            'makeblastdb': dict(VERSIONS['makeblastdb'], version='fixture-2')})
        self.assertEqual(len({first, changed, upgraded}), 3)
        self.assertTrue(first.with_suffix('.psq').exists())

    def test_corruption_is_rejected_and_never_overwritten(self):
        db, _ = self.prepare(FixtureRunner())
        broken = db.with_suffix('.psq')
        broken.write_bytes(b'corrupt')
        runner = FixtureRunner()
        with self.assertRaisesRegex(ValueError, 'Cache evidence changed'):
            self.prepare(runner)
        self.assertEqual(broken.read_bytes(), b'corrupt')
        self.assertEqual(runner.commands, [])

    def test_interrupted_build_is_not_published_and_can_be_retried(self):
        with self.assertRaisesRegex(RuntimeError, 'Interrupted'):
            self.prepare(FixtureRunner(fail=True))
        self.assertEqual(list(self.cache.rglob('cache.json')), [])
        self.assertEqual(list(self.cache.rglob('*.building-*')), [])
        _, receipt = self.prepare(FixtureRunner())
        self.assertFalse(receipt['reused'])

    def test_unsealed_alias_is_rejected_before_search(self):
        db, _ = self.prepare(FixtureRunner())
        db.with_suffix('.pal').write_text('DBLIST unverified-database\n')
        with self.assertRaisesRegex(ValueError, 'Unexpected shared search-cache'):
            self.prepare(FixtureRunner())

    def test_simultaneous_queries_publish_only_one_database(self):
        runners = [FixtureRunner() for _ in range(4)]
        with ThreadPoolExecutor(max_workers=4) as pool:
            values = list(pool.map(self.prepare, runners))
        self.assertEqual(len({path for path, _ in values}), 1)
        self.assertEqual(sum(len(r.commands) for r in runners), 1)
        self.assertEqual(sum(not r['reused'] for _, r in values), 1)

    def test_source_change_during_build_does_not_publish_stale_key(self):
        original = self.fasta.read_text()
        runner = FixtureRunner()
        original_run = runner.run

        def mutate(*args):
            original_run(*args)
            self.fasta.write_text(original+'A\n')

        runner.run = mutate
        with self.assertRaisesRegex(ValueError, 'FASTA changed'):
            self.prepare(runner)
        self.assertEqual(list(self.cache.rglob('cache.json')), [])


@unittest.skipUnless(shutil.which('samtools'), 'requires actual samtools')
class GenomeCache(unittest.TestCase):
    def test_gzip_index_is_reused_and_owned_fasta_survives_source_edits(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            source = root/'genome.fa.gz'
            with gzip.open(source, 'wt') as out:
                out.write('>chr\nACGTacgtNNNN\n')
            cold = Runner(root/'cold')
            genome, receipt = prepare_genome(source, root/'cache', cold)
            self.assertEqual(genome.fetch('chr', 0, 12, preserve_case=True), 'ACGTacgtNNNN')
            owned = genome.path
            genome.close()
            self.assertFalse(receipt['reused'])
            warm = Runner(root/'warm')
            genome, receipt = prepare_genome(source, root/'cache', warm)
            genome.close()
            self.assertTrue(receipt['reused'])
            self.assertEqual(warm.commands, [])
            self.assertEqual(len(cold.commands), 1)
            self.assertTrue(Path(cold.commands[0]['stderr']).exists())
            with gzip.open(source, 'wt') as out:
                out.write('>chr\nTTTTacgtNNNN\n')
            new, receipt = prepare_genome(source, root/'cache', Runner(root/'new'))
            try:
                self.assertNotEqual(new.path, owned)
                self.assertEqual(owned.read_text(), '>chr\nACGTacgtNNNN\n')
            finally:
                new.close()


@unittest.skipUnless(all(shutil.which(x) for x in ['samtools', 'makeblastdb', 'blastp', 'tblastn']),
                     'requires actual BLAST and samtools')
class SearchEquivalence(unittest.TestCase):
    def test_annotation_and_genomic_routes_keep_identical_cold_warm_results(self):
        from genecontext.cli import parser
        from genecontext.sequence import CODE, write_fasta, write_tsv
        from genecontext.stage1 import search, verify
        rng = random.Random(814)
        protein = 'M'+''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(99))
        codons = {aa: next(c for c, x in CODE.items() if x == aa) for aa in protein}
        cds = ''.join(codons[aa] for aa in protein)+'TAA'
        flank = lambda: ''.join(rng.choice('ACGT') for _ in range(600))
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            write_fasta(root/'genome.fna', [('chr', flank()+cds+flank())])
            write_fasta(root/'query.faa', [('reference', protein)])
            write_fasta(root/'proteins.faa', [('XP_1', protein)])
            (root/'genes.gff3').write_text(
                '##gff-version 3\n##sequence-region chr 1 1503\n'
                'chr\tx\tgene\t601\t903\t.\t+\t.\tID=g;Name=TEST\n'
                'chr\tx\tmRNA\t601\t903\t.\t+\t.\tID=tx;Parent=g\n'
                'chr\tx\tCDS\t601\t903\t.\t+\t0\tParent=tx;protein_id=XP_1\n')
            for route in ['genome', 'annotated']:
                with self.subTest(route=route):
                    row = dict(genome_id='species', fasta=str(root/'genome.fna'))
                    if route == 'annotated':
                        row.update(gff3=str(root/'genes.gff3'), protein=str(root/'proteins.faa'))
                    manifest = root/(route+'.tsv')
                    write_tsv(manifest, [row], list(row))
                    common = ['stage1', '--gene', 'TEST', '--query', str(root/'query.faa'),
                              '--query-type', 'protein', '--genomes', str(manifest),
                              '--neighbor-cache', str(root/(route+'_cache')), '--early-synteny', 'off',
                              '--refine', 'none', '--threads', '1']
                    values = []
                    for label in ['cold', 'warm']:
                        output = root/(route+'_'+label)
                        self.assertEqual(search(parser().parse_args(common+['--out', str(output)])), 0)
                        values.append(verify(output)[1])
                    self.assertTrue(values[0]['loci'])
                    self.assertEqual(values[0]['loci'], values[1]['loci'])
                    self.assertFalse(values[0]['genomes'][0]['search_database_cache']['reused'])
                    self.assertTrue(values[1]['genomes'][0]['search_database_cache']['reused'])
                    self.assertEqual(len(values[0]['queries']), len(values[1]['queries']))
                    tools = [Path(c['argv'][0]).name for c in values[1]['commands']]
                    self.assertNotIn('makeblastdb', tools)
                    self.assertNotIn('samtools', tools)
                    for suffix in ['blast.tsv', 'hsps.json']:
                        self.assertEqual((root/(route+'_cold')/'raw'/('species.'+suffix)).read_bytes(),
                                         (root/(route+'_warm')/'raw'/('species.'+suffix)).read_bytes())
