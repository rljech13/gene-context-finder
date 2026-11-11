"""Bounded software stress controls; never biological review decisions."""
import gzip
import json
from pathlib import Path
import random
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.sequence import Model, coding_sequence, sha256, write_fasta, write_tsv
from genecontext.stage2 import extract_intervals, neighbor_pairs, parse_hits, run, selected_regions
from genecontext.stage2_db import prepare, verify_db


class BoundedStress(unittest.TestCase):
    def test_nonboolean_commit_cannot_approve_draft_selections(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'draft.json'
            for value in ['false',1,None]:
                path.write_text(json.dumps(dict(schema='stage1_pair_draft_v1',review_id='x',pairs=[dict(genome_id='g',committed=value)])))
                args=SimpleNamespace(selection=None,draft=str(path),committed_only=True)
                with self.subTest(value=value),self.assertRaisesRegex(ValueError,'boolean'):
                    selected_regions(dict(review_id='x',genomes=[dict(genome_id='g')]),args)

    def test_8192_overlapping_intervals_from_16_mib_gzip(self):
        pattern='ACGTrymkSWBDHVN';size=4*1024*1024
        dna=(pattern*((size+len(pattern)-1)//len(pattern)))[:size]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'assembly.fa.gz'
            with gzip.open(path,'wt') as out:
                for c in range(4):
                    out.write('>chr'+str(c)+'\n')
                    for i in range(0,size,127):out.write(dna[i:i+127]+'\n')
            requests={str(i):('chr'+str(i%4),(i*7919)%(size-129),(i*7919)%(size-129)+129) for i in range(8192)}
            result,lengths=extract_intervals(path,requests)
            self.assertEqual(len(result),8192)
            self.assertEqual(lengths,{'chr'+str(i):size for i in range(4)})
            for key,(_,a,b) in requests.items():self.assertEqual(result[key],dna[a:b])

    def test_76800_isoform_hits_collapse_into_gene_pairs_and_keep_competition(self):
        models=[dict(model_id=f'{g}_{k}_{i}',genome_id=str(g),gene_group=f'{g}_{k}',overlaps_anchor=g==0 and k==0)
                for g in range(4) for k in range(2) for i in range(40)]
        hits=[dict(query=a['model_id'],subject=b['model_id'],bitscore=100,evalue=1e-20,
                   identity=.8,coverage_query=1,coverage_subject=1)
              for a in models for b in models if a['genome_id']!=b['genome_id']]
        self.assertEqual(len(hits),76800)
        args=SimpleNamespace(identification_coverage=.5,identification_identity=.3,competing_ratio=.95)
        pairs=neighbor_pairs(models,hits,args)
        self.assertEqual(len(pairs),48)
        self.assertTrue(all(p['status']=='ambiguous_proposal' for p in pairs))
        self.assertEqual(sum(p['kind']=='target_similarity' for p in pairs),12)
        expected={(p['gene_a'],p['gene_b'],p['status'],p['kind']) for p in pairs}
        random.Random(7).shuffle(hits)
        self.assertEqual(expected,{(p['gene_a'],p['gene_b'],p['status'],p['kind']) for p in neighbor_pairs(models,hits,args)})

    def test_256_random_spliced_cds_both_strands_and_separate_stops(self):
        rng=random.Random(19);codons={'A':'GCT','D':'GAT','Q':'CAG','G':'GGT'}
        def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
        for i in range(256):
            aa='M'+''.join(rng.choice(list(codons)) for _ in range(rng.randrange(5,60)))
            cds='ATG'+''.join(codons[a] for a in aa[1:])+'TAA'
            cut1=rng.randrange(1,len(cds)-7);cut2=rng.randrange(cut1+1,len(cds)-3)
            slices=[cds[:cut1],cds[cut1:cut2],cds[cut2:]];gene='';blocks=[];used=0
            for part in slices:
                if gene:gene+='G'*rng.randrange(4,29)
                blocks.append((len(gene),len(gene)+len(part),(3-used%3)%3));gene+=part;used+=len(part)
            strand='-' if i%2 else '+';span=len(gene)
            if strand=='-':gene=rc(gene);blocks=[(span-b,span-a,p) for a,b,p in blocks]
            stops=[]
            if i%3==0:
                a,b,p=blocks[-1]
                if strand=='+':blocks[-1]=(a,b-3,p);stops=[(b-3,b)]
                else:blocks[-1]=(a+3,b,p);stops=[(a,a+3)]
            class Genome:
                def fetch(self,c,a,b,s='+'):return gene[a:b] if s=='+' else rc(gene[a:b])
            model=Model(str(i),'chr',strand,cds=blocks,stops=stops,gene='g')
            observed,protein,status,issues=coding_sequence(model,Genome())
            self.assertEqual((observed,protein,status,issues),(cds,aa,'complete_predicted',[]))

    def test_blast_rejects_nonfinite_or_out_of_range_statistics(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'hits.tsv'
            for index,value in [(2,'nan'),(2,'101'),(2,'-1'),(8,'nan'),(8,'-1'),(9,'inf'),(9,'-1')]:
                row=['q','s','50','60','1','60','1','60','1e-20','100','100','100'];row[index]=value
                path.write_text('\t'.join(row)+'\n')
                with self.subTest(column=index,value=value),self.assertRaises(ValueError):parse_hits(path,{'q':'A'*100},{'s':100})


@unittest.skipUnless(shutil.which('blastp') and shutil.which('makeblastdb'),'requires BLAST+')
class FailureAndIntegrity(unittest.TestCase):
    def test_reference_corruption_and_extra_blast_volume_are_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);write_fasta(root/'refs.fa',[('ref','MACDEFGHIKLMNPQRSTVWY'*4)])
            prepare(SimpleNamespace(fasta=str(root/'refs.fa'),out=str(root/'db')))
            verify_db(root/'db');protein=root/'db/proteins.faa';old=protein.read_text();protein.write_text(old+'A')
            with self.assertRaisesRegex(ValueError,'changed'):verify_db(root/'db')
            protein.write_text(old);(root/'db/blastdb.unexpected').write_text('bad')
            with self.assertRaisesRegex(ValueError,'file set changed'):verify_db(root/'db')

    def test_prediction_failure_preserves_annotation_and_reports_failed_region(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);protein='M'+'ACDEFGHIKLMNPQRSTVWY'*6
            # A long complete CDS whose identity is immaterial to the deliberate
            # predictor-failure control. The BLAST tools themselves are real.
            cds='ATG'+'GCT'*120+'TAA';dna='C'*1000+cds+'C'*1000
            write_fasta(root/'genome.fa',[('chr',dna)]);write_fasta(root/'refs.fa',[('ref',protein)])
            prepare(SimpleNamespace(fasta=str(root/'refs.fa'),out=str(root/'db')))
            (root/'genes.gff').write_text('##gff-version 3\n'+
                f'chr\ttest\tgene\t1001\t{1000+len(cds)}\t.\t+\t.\tID=g;Name=TEST\n'+
                f'chr\ttest\tmRNA\t1001\t{1000+len(cds)}\t.\t+\t.\tID=t;Parent=g\n'+
                f'chr\ttest\tCDS\t1001\t{1000+len(cds)}\t.\t+\t0\tParent=t\n')
            source=root/'stage1';source.mkdir();data=dict(review_id='stress',config=dict(gene='TEST',window=500,query_type='protein'),
                queries=[dict(id='q',sequence=protein,members=[dict(query='control')])],genomes=[dict(genome_id='a',status='success',fasta_sha256=sha256(root/'genome.fa'))],
                loci=[dict(genome_id='a',locus_id='l',gene='TEST',contig='chr',start=1000,end=1000+len(cds),strand='+',contig_length=len(dna),pass_filters=True,boundary_source='control')])
            (source/'evidence.json').write_text(json.dumps(data))
            write_tsv(root/'datasets.tsv',[dict(genome_id='a',fasta=str(root/'genome.fa'),gff3=str(root/'genes.gff'))],['genome_id','fasta','gff3'])
            write_tsv(root/'decisions.tsv',[dict(review_id='stress',genome_id='a',locus_id='l',decision='inspect',reason='check_context',notes='software control only')],['review_id','genome_id','locus_id','decision','reason','notes'])
            args=parser().parse_args(['stage2','--results',str(source),'--decisions',str(root/'decisions.tsv'),'--genomes',str(root/'datasets.tsv'),'--protein-db',str(root/'db'),'--out',str(root/'result')])
            def failed(args,out,rows,windows,runner,errors):
                errors.append(dict(stage='context_prediction',locus_id='l',error='Deliberate software failure control'))
                return {},dict(batches=[dict(loci=['l'],status='error')])
            with patch('genecontext.stage2.verify',return_value=(source,data)),patch('genecontext.stage2.validate_prediction_inputs',return_value=None),patch('genecontext.tiberius_backend.predict_batches',side_effect=failed):
                self.assertEqual(run(args),1)
            result=json.loads((root/'result/evidence.json').read_text())
            self.assertEqual(result['status'],'completed_with_errors');self.assertEqual(len(result['models']),1)
            self.assertEqual(result['models'][0]['source'],'annotation');self.assertTrue((root/'result/report.html').is_file())
            self.assertTrue(result['errors']);self.assertTrue((root/'result/windows.fna').is_file())
            # A truncated ordinary gzip is an extraction failure, not an
            # uncaught EOFError or a successful empty prediction.
            with gzip.open(root/'truncated.fa.gz','wt') as out:out.write('>chr\n'+dna+'\n')
            raw=(root/'truncated.fa.gz').read_bytes();(root/'truncated.fa.gz').write_bytes(raw[:-8])
            data['genomes'][0].pop('fasta_sha256')
            write_tsv(root/'datasets.tsv',[dict(genome_id='a',fasta=str(root/'truncated.fa.gz'),gff3=str(root/'genes.gff'))],['genome_id','fasta','gff3'])
            args.out=str(root/'truncated_result')
            with patch('genecontext.stage2.verify',return_value=(source,data)),patch('genecontext.stage2.validate_prediction_inputs',return_value=None),patch('genecontext.tiberius_backend.predict_batches') as predictor:
                self.assertEqual(run(args),1);predictor.assert_not_called()
            bad=json.loads((root/'truncated_result/evidence.json').read_text())
            self.assertEqual(bad['status'],'completed_with_errors');self.assertEqual(bad['models'],[])
            self.assertEqual(bad['provenance']['prediction']['status'],'no_extracted_windows')
            self.assertTrue((root/'truncated_result/report.html').is_file())
