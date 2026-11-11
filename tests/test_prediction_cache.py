"""Prediction reuse never substitutes a different DNA/window/method or failed run."""
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.sequence import read_fasta, write_fasta
from genecontext.tiberius_backend import predict_batches
from genecontext.tools import Runner


class FixtureRunner(Runner):
    broken = False
    def run(self, args, output, log=None):
        self.commands.append(dict(argv=list(map(str,args)),output=str(output),stderr=str(log),returncode=0))
        fasta=Path(args[args.index('--fasta')+1]); dest=Path(args[args.index('--out')+1])
        locus=next(iter(read_fasta(fasta)))
        dest.write_text(f'##gff-version 3\n{locus}\tTiberius\tmRNA\t1\t303\t.\t+\t.\tID=t1;Parent=g1\n{locus}\tTiberius\tCDS\t1\t303\t.\t+\t0\tParent=t1\n')
        Path(str(dest)+'.raw.gff3').write_text(dest.read_text())
        write_fasta(str(dest)+'.cds.fna',[('x|t1|coding','ATG'+'GCT'*100)])
        write_fasta(str(dest)+'.proteins.faa',[('x|t1|protein', 'M'+'A'*100 if not self.broken else 'WRONG')])
        Path(output).write_text('GENE_CONTEXT_CPU='+json.dumps(dict(visible_devices=['CPU:0'],input_windows=1))+'\n')
        Path(log).write_text('')


class PredictionCache(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.args=SimpleNamespace(tiberius_python='fixture',predictor_cache=str(self.root/'cache'),
            offline=True,threads=1,tiberius_batch_size=1,tiberius_worker_mode='isolated')
        self.runtime=dict(packages=dict(tiberius='2.0.9'),python='control')
        self.weights=dict(path='control',files_sha256=dict(weights='verified-control'))
        self.probe=patch('genecontext.tiberius_backend.runtime_info',return_value=self.runtime);self.probe.start()
        self.model=patch('genecontext.tiberius_backend.ensure_model',return_value=self.weights);self.model.start()
        self.addCleanup(self.probe.stop);self.addCleanup(self.model.stop);self.addCleanup(self.temp.cleanup)
        self.counter=0

    def predict(self, locus='window', sequence=None, broken=False):
        self.counter+=1;root=self.root/str(self.counter);root.mkdir();runner=FixtureRunner(root);runner.broken=broken
        sequence=sequence or 'ATG'+'GCT'*100+'N'*1000
        rows=[dict(locus_id=locus,context_predictor_model='vertebrates',window_bp=len(sequence))];errors=[]
        result,proof=predict_batches(self.args,root,rows,{locus:sequence},runner,errors)
        return result,proof,runner,errors

    def test_repeated_window_skips_worker_and_roundtrips_identically(self):
        first,receipt,runner,errors=self.predict();self.assertEqual(errors,[]);self.assertEqual(len(runner.commands),1)
        second,receipt,runner,errors=self.predict();self.assertEqual(errors,[]);self.assertEqual(first,second)
        self.assertEqual(runner.commands,[]);self.assertTrue(receipt['batches'][0]['cache']['reused'])

    def test_each_input_and_method_change_invalidates_cache(self):
        self.predict()
        for change in ['dna','header','threads','batch','weights','runtime']:
            kwargs={}
            if change=='dna':kwargs['sequence']='ATG'+'GCT'*100+'A'+'N'*999
            if change=='header':kwargs['locus']='other'
            if change=='threads':self.args.threads=2
            if change=='batch':self.args.tiberius_batch_size=2
            if change=='weights':self.weights['files_sha256']['weights']='changed'
            if change=='runtime':self.runtime['python']='changed'
            _,proof,runner,errors=self.predict(**kwargs)
            self.assertEqual(errors,[],change);self.assertEqual(len(runner.commands),1,change)
            self.assertFalse(proof['batches'][0]['cache']['reused'],change)

    def test_invalid_roundtrip_is_never_published(self):
        result,proof,runner,errors=self.predict(broken=True)
        self.assertEqual(result,{});self.assertIn('roundtrip',errors[0]['error'])
        self.assertFalse(list((self.root/'cache/predictions').glob('*/cache.json')))
        _,proof,runner,errors=self.predict();self.assertEqual(errors,[]);self.assertEqual(len(runner.commands),1)

    def test_tampered_cache_is_preserved_and_explicitly_rejected(self):
        _,proof,_,_=self.predict();stored=Path(proof['batches'][0]['cache']['path'])/'predictions.gff3'
        stored.write_text('tampered')
        result,proof,runner,errors=self.predict();self.assertEqual(result,{})
        self.assertIn('Cache evidence changed',errors[0]['error']);self.assertFalse(runner.commands)
        self.assertEqual(stored.read_text(),'tampered')

    def test_valid_empty_prediction_is_cacheable(self):
        original=FixtureRunner.run
        def empty(runner,args,output,log=None):
            original(runner,args,output,log);dest=Path(args[args.index('--out')+1])
            for suffix in ['', '.raw.gff3','.cds.fna','.proteins.faa']:Path(str(dest)+suffix).write_text('')
        with patch.object(FixtureRunner,'run',empty):
            result,_,_,errors=self.predict();self.assertEqual(errors,[]);self.assertEqual(result,{'window':[]})
            result,proof,runner,errors=self.predict();self.assertEqual(errors,[]);self.assertFalse(runner.commands)
