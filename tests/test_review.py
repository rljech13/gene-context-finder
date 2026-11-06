"""Review gates, immutable evidence, and actual staged search/prediction/comparison."""
import csv
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.comparison import table
from genecontext.review import contexts, decisions, review_search, verify, validate_model_sequence
from genecontext.sequence import read_fasta, reverse_complement, sha256, translate, write_fasta, write_tsv


def edit(path, change):
    rows=table(path)
    for row in rows:change(row)
    write_tsv(path,rows,list(rows[0]))


class ReviewUnits(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.source=self.root/'search';self.source.mkdir()
        queries=[dict(id='q0001',name='opaque',protein='M'+'A'*99,coding='',type='protein')]
        original=dict(queries=queries,genomes=[dict(genome_id='g'),dict(genome_id='missing')],status='completed',errors=[])
        (self.source/'manifest.json').write_text(json.dumps(original))
        base=dict(query='q0001',genome_id='g',contig='c',strand='+',status='complete_predicted',source='annotation',protein_aa=100,coverage=1,identity=1,cds_bp=303)
        rows=[dict(base,model_id='first',source_model='a.t',source_gene='a',start=100,end=403,cds_exons='[[100,403,0]]'),dict(base,model_id='second',source_model='b.t',source_gene='b',start=1000,end=1303,cds_exons='[[1000,1303,0]]')]
        write_tsv(self.source/'models.tsv',rows,list(rows[0]));write_tsv(self.source/'dna_regions.tsv',[],['region_id','query','genome_id'])
        seq=self.source/'sequences/q0001';seq.mkdir(parents=True)
        write_fasta(seq/'protein.faa',[(r['model_id'],'M'+'A'*99) for r in rows]);write_fasta(seq/'cds.fna',[(r['model_id'],'ATG'+'GCT'*99+'TAA') for r in rows])
        self.review=self.root/'review';review_search(parser().parse_args(['review-search','--results',str(self.source),'--out',str(self.review)]))

    def tearDown(self):self.tmp.cleanup()

    def args(self):
        return parser().parse_args(['contexts','--review',str(self.review),'--decisions',str(self.review/'loci.decisions.tsv'),'--out',str(self.root/'contexts')])

    def choose(self):
        edit(self.review/'loci.decisions.tsv',lambda r:r.update(decision='inspect' if r['genome_id']=='g' else 'unresolved',model_ids='first;second' if r['genome_id']=='g' else '',notes='No candidate in the control' if r['genome_id']=='missing' else 'Inspect both copies'))

    def test_blank_decisions_block_before_prediction_or_output_creation(self):
        with patch('genecontext.comparison.compare',side_effect=AssertionError('must not execute')):
            with self.assertRaisesRegex(ValueError,'Pending'):contexts(self.args())
        self.assertFalse((self.root/'contexts').exists())

    def test_multiple_explicit_loci_advance_and_unresolved_pair_is_retained(self):
        self.choose()
        with patch('genecontext.comparison.compare',return_value=0) as call:
            self.assertEqual(contexts(self.args()),0)
        selected,_,observations=call.call_args.kwargs['selection']
        self.assertEqual([r['model_id'] for r in selected],['first','second'])
        self.assertEqual(len({r['locus_id'] for r in selected}),2)
        self.assertEqual(observations[1]['status'],'review_unresolved')
        self.assertTrue(call.call_args.kwargs['context_only'])

    def test_wrong_pair_ids_duplicate_rows_and_missing_decisions_are_rejected(self):
        self.choose();edit(self.review/'loci.decisions.tsv',lambda r:r.update(model_ids='first') if r['genome_id']=='missing' else None)
        with self.assertRaisesRegex(ValueError,'Only inspect'):contexts(self.args())
        self.choose();p=self.review/'loci.decisions.tsv';rows=table(p);write_tsv(p,[rows[0],rows[0]],list(rows[0]))
        with self.assertRaisesRegex(ValueError,'duplicate'):contexts(self.args())
        write_tsv(p,[rows[0]],list(rows[0]))
        with self.assertRaisesRegex(ValueError,'Missing explicit'):contexts(self.args())

    def test_modifying_evidence_is_detected_but_decision_edits_are_allowed(self):
        self.choose();verify(self.review,'locus_review')
        p=self.review/'candidates.tsv';p.write_text(p.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'evidence changed'):verify(self.review,'locus_review')

    def test_decisions_from_another_checkpoint_cannot_be_replayed(self):
        self.choose();other=self.root/'another_review'
        review_search(parser().parse_args(['review-search','--results',str(self.source),'--out',str(other)]))
        args=self.args();args.review=str(other)
        with self.assertRaisesRegex(ValueError,'another review checkpoint'):contexts(args)

    def test_relocated_checkpoint_and_source_keep_the_review_identity(self):
        self.choose();moved=self.root/'moved_review';shutil.copytree(self.review,moved)
        relocated=self.root/'moved_search';self.source.rename(relocated)
        args=self.args();args.review=str(moved);args.decisions=str(moved/'loci.decisions.tsv');args.search_results=str(relocated)
        with patch('genecontext.comparison.compare',return_value=0) as call:
            self.assertEqual(contexts(args),0)
        self.assertEqual(call.call_args.args[0].results,str(relocated.resolve()))
        self.assertEqual(json.loads((moved/'manifest.json').read_text())['review_id'],json.loads((self.review/'manifest.json').read_text())['review_id'])

    def test_source_export_changes_are_detected(self):
        p=self.source/'sequences/q0001/protein.faa';p.write_text(p.read_text().replace('AAA','CCC'))
        with self.assertRaisesRegex(ValueError,'Source search evidence changed'):verify(self.review,'locus_review')

    def test_explicitly_rejecting_all_pairs_finishes_without_computation(self):
        edit(self.review/'loci.decisions.tsv',lambda r:r.update(decision='unresolved',notes='Deliberately unresolved control'))
        with patch('genecontext.comparison.compare',side_effect=AssertionError('must not execute')):
            self.assertEqual(contexts(self.args()),0)
        manifest=json.loads((self.root/'contexts/manifest.json').read_text())
        self.assertEqual(manifest['status'],'completed_no_comparable_targets')
        self.assertEqual(manifest['commands'],[])

    def test_chosen_model_requires_native_dna_roundtrip_on_both_strands(self):
        cds='ATG'+'GCT'*20+'TAA';window='C'*100+cds+'C'*100
        for strand,seq in [('+',window),('-',window)]:
            row=dict(strand=strand,contig='c',window_start=1000,window_end=1000+len(window))
            a,b=(1100,1100+len(cds)) if strand=='+' else (row['window_end']-100-len(cds),row['window_end']-100)
            option=dict(option_id='chosen',model_id='m',strand=strand,native_cds=[[a,b,0]])
            validate_model_sequence(option,row,seq,dict(cds=cds,protein=translate(cds).removesuffix('*')))
            with self.assertRaisesRegex(ValueError,'round-trip'):
                validate_model_sequence(option,row,seq,dict(cds=cds[:-3]+'TAG',protein=translate(cds).removesuffix('*')))


TOOLS=['samtools','miniprot','tblastn','blastn','blastp','makeblastdb','mafft']
@unittest.skipUnless(all(shutil.which(x) for x in TOOLS),'requires actual search/alignment tools')
class StagedTiberiusIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from genecontext.tiberius_backend import runtime_info,runtime_python,ensure_model
        args=parser().parse_args(['compare','--results','unused','--out','unused'])
        cls.python=os.environ.get('GENE_CONTEXT_TIBERIUS_PYTHON') or runtime_python(args)
        try:runtime_info(cls.python);ensure_model('vertebrates',args.predictor_cache,offline=True)
        except (ValueError,OSError) as e:raise unittest.SkipTest(str(e))

    def test_manual_stops_prediction_model_selection_and_predictor_free_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures=Path(__file__).parent/'fixtures';native=next(iter(read_fasta(fixtures/'HS04636.fna').values()));native=native+'N'*5000+native
            genome_rows=[]
            for gid,seq in [('forward',native),('reverse',reverse_complement(native))]:
                path=root/(gid+'.fna');write_fasta(path,[('chr',seq)]);genome_rows.append(dict(genome_id=gid,fasta=str(path)))
            inputs=root/'genomes.tsv';write_tsv(inputs,genome_rows,['genome_id','fasta'])
            search=root/'01';stage2=root/'02';out=root/'03'
            def run(command,expected=0):
                result=subprocess.run([sys.executable,'-m','genecontext',*map(str,command)],capture_output=True,text=True)
                self.assertEqual(result.returncode,expected,result.stdout+result.stderr)
                return result
            run(['search','--query',fixtures/'PTGS2.reference.faa','--query-type','protein','--genomes',inputs,'--out',search,'--threads','1'])
            gate=search/'review';self.assertEqual(json.loads((gate/'manifest.json').read_text())['status'],'awaiting_review')
            self.assertFalse((search/'comparison').exists())
            run(['contexts','--review',gate,'--decisions',gate/'loci.decisions.tsv','--out',stage2],2)
            self.assertFalse(stage2.exists())
            edit(gate/'loci.decisions.tsv',lambda r:r.update(decision='inspect',model_ids=r['suggested_model_id'],notes='Computational PTGS2 control'))
            run(['contexts','--review',gate,'--decisions',gate/'loci.decisions.tsv','--out',stage2,'--tiberius-python',self.python,'--offline','--context-reference','none','--threads','1'])
            self.assertFalse((stage2/'dna_pairs.tsv').exists());self.assertFalse((stage2/'queries').exists())
            state=json.loads((stage2/'state.json').read_text());proof=json.loads((stage2/'context_prediction.json').read_text())
            self.assertEqual(len(proof['batches']),2)
            options={}
            for r in state['observations']:
                matching=[o for o in state['options'] if o['genome_id']==r['genome_id'] and o['origin']=='prediction' and o['eligible_target'] and o['protein_aa']>=500]
                self.assertTrue(matching);seed=next(x for x in state['selected'] if x['genome_id']==r['genome_id']);options[r['genome_id']]=min(matching,key=lambda o:abs(o['start']-seed['target_start']))
            run(['compare-reviewed','--review',stage2,'--models',stage2/'models.decisions.tsv','--neighbors',stage2/'neighbors.decisions.tsv','--out',out],2)
            self.assertFalse(out.exists())
            edit(stage2/'models.decisions.tsv',lambda r:r.update(decision='keep',option_id=options[r['genome_id']]['option_id'],target_name='PTGS2_control',allow_partial='yes',allow_short_flanks='yes',notes='Explicit control: fixture ends truncate the requested flanks'))
            # A choice of targets alone cannot silently approve the other genes.
            neighbors=table(stage2/'neighbors.decisions.tsv');self.assertGreaterEqual(len(neighbors),4)
            run(['compare-reviewed','--review',stage2,'--models',stage2/'models.decisions.tsv','--neighbors',stage2/'neighbors.decisions.tsv','--out',out],2)
            self.assertFalse(out.exists())
            edit(stage2/'neighbors.decisions.tsv',lambda r:r.update(decision='unnamed',notes='Keep unnamed control predictions'))
            before={str(p.relative_to(stage2)):sha256(p) for p in stage2.rglob('predictions.gff3')}
            run(['compare-reviewed','--review',stage2,'--models',stage2/'models.decisions.tsv','--neighbors',stage2/'neighbors.decisions.tsv','--out',out,'--threads','1'])
            manifest=json.loads((out/'manifest.json').read_text());self.assertEqual(manifest['status'],'completed',manifest['errors'])
            self.assertTrue(manifest['review']['prediction_reused'])
            self.assertEqual({str(p.relative_to(stage2)):sha256(p) for p in stage2.rglob('predictions.gff3')},before)
            self.assertTrue(all(Path(c['argv'][0]).name in {'blastp','blastn','makeblastdb','mafft'} for c in manifest['commands']))
            aa=read_fasta(out/'reviewed_target_proteins.faa');cds=read_fasta(out/'reviewed_target_cds.fna')
            self.assertEqual(len(aa),2)
            for lid in aa:self.assertEqual(translate(cds[lid]).removesuffix('*'),aa[lid]);self.assertGreaterEqual(len(aa[lid]),500)
            self.assertTrue(all(r['review_origin']=='prediction' for r in table(out/'reviewed_targets.tsv')))
            self.assertEqual(len(read_fasta(out/'queries/q0001/proteins.aligned.faa')),3)
            self.assertEqual(len(table(out/'dna_pairs.tsv')),1)
            self.assertTrue((out/'approved_models.tsv').exists());self.assertTrue((out/'reviewed_targets.gff3').exists())


@unittest.skipUnless(all(shutil.which(x) for x in TOOLS),'requires actual search/alignment tools')
class StagedDNAIntegration(unittest.TestCase):
    def test_noncoding_review_without_any_coding_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rng=random.Random(551);dna=''.join(rng.choices('ACGT',k=2000))
            query=root/'query.fna';write_fasta(query,[('noncoding',dna[500:901])])
            inputs=root/'genomes.tsv';rows=[]
            for gid in ['one','two']:
                path=root/(gid+'.fa');write_fasta(path,[('chr',dna)]);rows.append(dict(genome_id=gid,fasta=str(path)))
            write_tsv(inputs,rows,['genome_id','fasta'])
            def run(command):
                result=subprocess.run([sys.executable,'-m','genecontext',*map(str,command)],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            stage1=root/'01';stage2=root/'02';out=root/'03'
            run(['search','--query',query,'--query-type','dna','--genomes',inputs,'--out',stage1,'--threads','1'])
            gate=stage1/'review'
            edit(gate/'loci.decisions.tsv',lambda r:r.update(decision='inspect',model_ids=r['suggested_model_id']))
            run(['contexts','--review',gate,'--decisions',gate/'loci.decisions.tsv','--out',stage2,'--window','100','--context-mode','annotation','--context-reference','none','--threads','1'])
            options=table(stage2/'model_options.tsv');self.assertEqual({r['origin'] for r in options},{'dna_only'})
            lookup={r['genome_id']:r['option_id'] for r in options}
            edit(stage2/'models.decisions.tsv',lambda r:r.update(decision='keep',option_id=lookup[r['genome_id']]))
            run(['compare-reviewed','--review',stage2,'--models',stage2/'models.decisions.tsv','--neighbors',stage2/'neighbors.decisions.tsv','--out',out,'--threads','1'])
            self.assertEqual((out/'reviewed_target_proteins.faa').stat().st_size,0)
            self.assertEqual({r['status'] for r in table(out/'loci.tsv')},{'dna_region'})
            self.assertEqual(len(table(out/'dna_pairs.tsv')),1)
            gff=(out/'reviewed_targets.gff3').read_text();self.assertIn('\tregion\t',gff);self.assertNotIn('\tgene\t',gff)


if __name__=='__main__':unittest.main()
