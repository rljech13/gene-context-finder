import csv
import functools
import gzip
import hashlib
import http.server
import io
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.context_predict import prediction_model
from genecontext.sequence import read_fasta, read_gff, reverse_complement, translate, write_fasta
from genecontext.tiberius_backend import MODELS, ensure_model, runtime_info, runtime_python, predict_batches
from genecontext.tools import Runner


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


class TiberiusCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.files = {'model_config.json': b'{"softmasking":false}', 'model_layers.json': b'{}', 'weights.h5': b'weight checksum control'}
        self.make_archive()
        handler = functools.partial(QuietHandler, directory=str(self.root))
        self.server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True); self.thread.start()
        self.spec = dict(url=f'http://127.0.0.1:{self.server.server_port}/model.tar.gz',
                         md5=hashlib.md5(self.payload).hexdigest(), sha256=hashlib.sha256(self.payload).hexdigest(),
                         files={name: hashlib.sha256(data).hexdigest() for name, data in self.files.items()})

    def make_archive(self, traversal=False):
        with tarfile.open(self.root/'model.tar.gz', 'w:gz') as packed:
            for name, data in self.files.items():
                member = tarfile.TarInfo(('../' if traversal else 'weights/')+name)
                member.size = len(data); packed.addfile(member, io.BytesIO(data))
        self.payload = (self.root/'model.tar.gz').read_bytes()

    def tearDown(self):
        self.server.shutdown(); self.server.server_close(); self.thread.join(); self.temp.cleanup()

    def test_download_pinned_model_offline_reuse_and_corruption_preservation(self):
        with patch.dict(MODELS, fixture=self.spec):
            info = ensure_model('fixture', self.root/'cache')
            model = Path(info['path'])
            self.assertEqual((model/'weights.h5').read_bytes(), self.files['weights.h5'])
            self.assertEqual(ensure_model('fixture', self.root/'cache', offline=True)['files_sha256'], self.spec['files'])
            (model/'weights.h5').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'pinned file checksums'):
                ensure_model('fixture', self.root/'cache', offline=True)
            self.assertEqual((model/'weights.h5').read_bytes(), b'changed')

    def test_missing_offline_model_and_wrong_archive_hash_are_not_ready(self):
        with patch.dict(MODELS, fixture=dict(self.spec, sha256='0'*64)):
            with self.assertRaisesRegex(ValueError, 'not cached'):
                ensure_model('fixture', self.root/'cache', offline=True)
            with self.assertRaisesRegex(ValueError, 'pinned SHA-256'):
                ensure_model('fixture', self.root/'cache')
        self.assertFalse((self.root/'cache/models/fixture/model').exists())

    def test_archive_traversal_is_rejected(self):
        self.make_archive(traversal=True)
        spec = dict(self.spec, md5=hashlib.md5(self.payload).hexdigest(), sha256=hashlib.sha256(self.payload).hexdigest())
        with patch.dict(MODELS, fixture=spec):
            with self.assertRaisesRegex(ValueError, 'archive member'):
                ensure_model('fixture', self.root/'cache')
        self.assertFalse((self.root/'cache/models/fixture/model').exists())

    def test_model_map_and_venv_executable_are_preserved(self):
        model_dir = self.root/'local_model'; model_dir.mkdir()
        mapping = self.root/'models.tsv'; mapping.write_text('genome_id\tmodel\nfirst\tmammalia\nsecond\tlocal_model\n')
        args = parser().parse_args(['compare', '--results', 'unused', '--out', 'unused', '--prediction-models', str(mapping)])
        self.assertEqual(prediction_model(args, {'genome_id':'first'}), 'mammalia')
        self.assertEqual(prediction_model(args, {'genome_id':'second'}), str(model_dir.resolve()))
        mapping.write_text('genome_id\tmodel\nfirst\tmammalia\nfirst\tvertebrates\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            prediction_model(args, {'genome_id':'first'})
        binary = self.root/'venv_python'; binary.symlink_to(sys.executable)
        args.tiberius_python = str(binary)
        self.assertEqual(runtime_python(args), str(binary))


@unittest.skipUnless(all(shutil.which(p) for p in ['samtools','miniprot','blastn','blastp','blastx','makeblastdb','mafft']), 'requires actual sequence tools')
class TiberiusIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        args = parser().parse_args(['compare', '--results','unused','--out','unused'])
        cls.python = os.environ.get('GENE_CONTEXT_TIBERIUS_PYTHON') or runtime_python(args)
        try:
            runtime_info(cls.python)
            ensure_model('vertebrates', args.predictor_cache, offline=True)
            ensure_model('mammalia', args.predictor_cache, offline=True)
        except (ValueError, OSError) as error:
            raise unittest.SkipTest(str(error))

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        fixtures = Path(__file__).parent/'fixtures'
        human = next(iter(read_fasta(fixtures/'HS04636.fna').values()))
        chicken = next(iter(read_fasta(Path(__file__).parent.parent/'examples/queries/chicken_Pou5f3.cds.fna').values()))
        self.sequence = human + 'N'*5000 + chicken + 'N'*5000
        self.references = fixtures/'PTGS2.reference.faa'
        self.genomes = self.root/'genomes.tsv'
        records=[]
        for gid, seq in [('first',self.sequence),('second',reverse_complement(self.sequence))]:
            path=self.root/(gid+'.fna.gz')
            with gzip.open(path,'wt') as h:h.write('>chr\n'+seq+'\n')
            records.append((gid,str(path),''))
        fake=self.root/'annotation.gff3'
        fake.write_text('##gff-version 3\nchr\tfixture\tgene\t8000\t8500\t.\t+\t.\tID=FICTITIOUS;Name=FICTITIOUS\n')
        records[0]=(*records[0][:2],str(fake))
        self.genomes.write_text('genome_id\tfasta\tgff3\n'+'\n'.join('\t'.join(r) for r in records)+'\n')
        self.query=self.root/'query.faa'
        protein=next(iter(read_fasta(Path(__file__).parent.parent/'examples/queries/chicken_Pou5f3.protein.faa').values()))
        rng=random.Random(144)
        write_fasta(self.query,[('Pou5f3',protein),('negative_control','M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=179)))])

    def tearDown(self):
        self.temp.cleanup()

    def execute(self, command, expected=0):
        run=subprocess.run([sys.executable,'-m','genecontext',*map(str,command),
                            '--tiberius-python',self.python,'--offline','--threads','1'],capture_output=True,text=True)
        self.assertEqual(run.returncode,expected,run.stdout+run.stderr)

    def test_complete_protein_query_pipeline_with_isolated_windows_and_no_gff(self):
        out=self.root/'analysis'
        self.execute(['analyze','--query',self.query,'--query-type','protein','--genomes',self.genomes,'--out',out,
                      '--context-proteins',self.references,'--context-reference','none'])
        result=out/'comparison'; manifest=json.loads((result/'manifest.json').read_text())
        self.assertEqual(manifest['status'],'completed',manifest.get('errors'))
        self.assertEqual(len(table(result/'summary.tsv')),4)
        self.assertEqual(len(table(result/'loci.tsv')),2)
        prediction=manifest['context_prediction']
        self.assertEqual(prediction['predictor'],'tiberius')
        self.assertEqual(len(prediction['batches']),2)
        self.assertTrue(all(len(b['loci'])==1 for b in prediction['batches']))
        self.assertTrue(all(b['cpu_execution']['input_windows']==1 for b in prediction['batches']))
        self.assertTrue(all('GPU' not in d for b in prediction['batches'] for d in b['cpu_execution']['visible_devices']))
        neighbors=table(result/'neighbors.tsv')
        named=[r for r in neighbors if r['name']=='PTGS2']
        self.assertEqual({r['genome_id'] for r in named},{'first','second'})
        self.assertTrue(all(r['evidence']=='ab_initio_prediction' for r in neighbors))
        self.assertNotIn('FICTITIOUS',{r['name'] for r in neighbors})
        self.assertEqual(len(table(result/'provided_neighbors.tsv')),1)
        aligned=read_fasta(result/'queries/q0001/proteins.aligned.faa')
        self.assertEqual(len(aligned),3)
        self.assertTrue(all(float(r['identity'])==1 for r in table(result/'queries/q0001/protein_pairs.tsv')))
        self.assertTrue(all(float(r['identity'])==1 for r in table(result/'dna_pairs.tsv')))
        self.assertGreaterEqual(int(table(result/'synteny.tsv')[0]['unique_anchors']),1)
        # Native-coordinate exports must reproduce the extracted coding sequences on both strands.
        cds=read_fasta(result/'context_prediction_cds.fna'); proteins=read_fasta(result/'context_prediction_proteins.faa')
        for gid,native in [('first',self.sequence),('second',reverse_complement(self.sequence))]:
            models,_=read_gff(result/'context_predictions'/(gid+'.gff3'))
            for model in models:
                seq=''.join(native[a:b] if model.strand=='+' else reverse_complement(native[a:b]) for a,b,_ in model.cds)
                pid=model.id.removesuffix('.t1')
                self.assertEqual(seq.upper(),cds[pid])
                self.assertEqual(translate(seq.upper()).removesuffix('*'),proteins[pid])
        self.assertIn('Tiberius on CPU',(out/'report.html').read_text())

    def test_unrelated_long_window_does_not_change_existing_prediction(self):
        short = next(iter(read_fasta(Path(__file__).parent/'fixtures/HS04636.fna').values()))
        sequences = {'short_window': short, 'long_window': short+'N'*30000}
        def predict(ids, name):
            root=self.root/name;root.mkdir()
            args=parser().parse_args(['compare','--results','unused','--out',str(root),
                                     '--tiberius-python',self.python,'--offline','--threads','1'])
            rows=[dict(locus_id=key,context_predictor_model='vertebrates',window_bp=len(sequences[key])) for key in ids]
            errors=[]
            models,proof=predict_batches(args,root,rows,sequences,Runner(root),errors)
            self.assertEqual(errors,[])
            self.assertTrue(all(b['cpu_execution']['input_windows']==1 for b in proof['batches']))
            return [(m.contig,m.strand,m.cds,m.stops) for m in models['short_window']]
        alone=predict(['short_window'],'alone')
        self.assertTrue(alone)
        self.assertEqual(alone,predict(['long_window','short_window'],'with_longer_window'))

    def test_model_map_splits_batches_for_saved_dna_search(self):
        query=self.root/'seed.fna'; write_fasta(query,[('noncoding_seed',self.sequence[:400])])
        mapping=self.root/'models.tsv'; mapping.write_text('genome_id\tmodel\nfirst\tvertebrates\nsecond\tmammalia\n')
        out=self.root/'mixed'
        self.execute(['analyze','--query',query,'--query-type','dna','--genomes',self.genomes,'--out',out,
                      '--prediction-models',mapping,'--context-proteins',self.references,'--context-reference','none'])
        manifest=json.loads((out/'comparison/manifest.json').read_text())
        self.assertEqual({b['model'] for b in manifest['context_prediction']['batches']},{'vertebrates','mammalia'})
        self.assertTrue(all(b['status']=='completed' for b in manifest['context_prediction']['batches']))
        self.assertEqual({r['context_predictor_model'] for r in table(out/'comparison/loci.tsv')},{'vertebrates','mammalia'})

    def test_identification_error_keeps_unnamed_predictions_and_short_windows_are_explicit(self):
        query=self.root/'seed.fna'; write_fasta(query,[('noncoding_seed',self.sequence[:400])])
        bad=self.root/'bad.faa'; bad.write_text('>invalid\n'+'M'*40+'*'+'A'*40+'\n')
        out=self.root/'failed_identification'
        self.execute(['analyze','--query',query,'--query-type','dna','--genomes',self.genomes,'--out',out,
                      '--context-proteins',bad,'--context-reference','none'],expected=1)
        neighbors=table(out/'comparison/neighbors.tsv')
        self.assertTrue(neighbors)
        self.assertTrue(all(r['identification_status']=='identification_error' for r in neighbors))
        self.assertTrue((out/'comparison/dna_pairs.tsv').exists())
        self.assertTrue((out/'report.html').exists())
        short=self.root/'short'
        self.execute(['compare','--results',out/'search','--out',short,'--window','0','--context-reference','none'])
        self.assertEqual({r['context_status'] for r in table(short/'loci.tsv')},{'prediction_window_too_short'})
        self.assertEqual(json.loads((short/'manifest.json').read_text())['context_prediction']['batches'],[])


if __name__=='__main__':
    unittest.main()
