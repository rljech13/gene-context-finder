"""Regression: chicken NANOG must not vanish when only its domain aligns."""
import csv
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

from genecontext.comparison import rank_key
from genecontext.protein_discovery import hsp_summary
from genecontext.sequence import read_fasta, reverse_complement, write_fasta


def table(path):
    with open(path) as handle:return list(csv.DictReader(handle,delimiter='\t'))


class RescueUnits(unittest.TestCase):
    def test_overlapping_hsps_do_not_inflate_score_or_coverage(self):
        a=dict(qstart=0,qend=60,start=10,end=70,bitscore=100,aligned_aa=60,identity=0.5,evalue=1e-20)
        b=dict(a,qstart=40,qend=100,start=200,end=260,bitscore=90)
        c=dict(a,qstart=60,qend=100,start=70,end=110,bitscore=60,aligned_aa=40)
        result=hsp_summary([a,b,c],200)
        self.assertEqual(result['bitscore'],160)
        self.assertEqual(result['coverage'],0.5)
        self.assertEqual(result['aligned_aa'],100)

    def test_normalized_protein_score_precedes_forced_long_alignment(self):
        wrong=dict(coverage=0.95,identity=0.5,status='complete_predicted',protein_aa=300,model_id='wrong',ranking_score=0.1)
        correct=dict(wrong,coverage=0.22,identity=0.57,model_id='correct',ranking_score=0.2)
        self.assertLess(rank_key(correct),rank_key(wrong))
        self.assertLess(rank_key(dict(wrong,ranking_score='')),rank_key(dict(correct,ranking_score='')))


@unittest.skipUnless(all(shutil.which(x) for x in ['samtools','miniprot','blastp','tblastn','makeblastdb']), 'requires real sequence tools')
class NanogRegression(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        fixtures=Path(__file__).parent/'fixtures'
        self.query=self.root/'query.faa'
        protein=next(iter(read_fasta(fixtures/'chicken_nanog.faa').values()))
        rng=random.Random(733)
        write_fasta(self.query,[('opaque_reference',protein),('negative_control','M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=308)))])
        sequence=next(iter(read_fasta(fixtures/'Danio_nanog_locus.fna').values()));n=len(sequence)
        rows=[]
        for gid,rev in [('forward',False),('reverse',True)]:
            fasta=self.root/(gid+'.fna');write_fasta(fasta,[('chr',reverse_complement(sequence) if rev else sequence)])
            gff=self.root/(gid+'.gff3');lines=['##gff-version 3']
            for l in (fixtures/'Danio_nanog_locus.gff3').read_text().splitlines():
                if l.startswith('#'):continue
                f=l.split('\t');f[0]='chr'
                # Names cannot supply the answer: erase all NANOG strings.
                f[8]=f[8].replace('nanog','opaque').replace('NANOG','opaque')
                if rev:
                    a,b=int(f[3])-1,int(f[4]);f[3],f[4]=str(n-b+1),str(n-a);f[6]='-' if f[6]=='+' else '+'
                lines.append('\t'.join(f))
            gff.write_text('\n'.join(lines)+'\n');rows.append((gid,str(fasta),str(gff)))
        self.annotated=self.root/'annotated.tsv';self.annotated.write_text('genome_id\tfasta\tgff3\n'+'\n'.join('\t'.join(x) for x in rows)+'\n')
        self.unannotated=self.root/'unannotated.tsv';self.unannotated.write_text('genome_id\tfasta\tgff3\n'+'\n'.join('\t'.join(x[:2])+'\t' for x in rows)+'\n')

    def tearDown(self):self.temp.cleanup()

    def run_search(self,label,manifest,extra=()):
        out=self.root/label
        result=subprocess.run([sys.executable,'-m','genecontext','run','--query',str(self.query),'--query-type','protein',
            '--genomes',str(manifest),'--out',str(out),'--threads','1','--skip-alignments',*extra],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        return out

    def test_annotated_nanog_recovered_independently_on_both_strands(self):
        fast=self.run_search('old',self.annotated,['--protein-search','miniprot'])
        self.assertEqual(table(fast/'models.tsv'),[])
        out=self.run_search('rescued',self.annotated)
        models=table(out/'models.tsv');self.assertEqual(len(models),2,models)
        self.assertEqual({m['strand'] for m in models},{'+','-'})
        self.assertTrue(all(m['source']=='annotation' and m['status']=='complete_predicted' for m in models))
        self.assertTrue(all(m['homology_support']=='local_alignment' and float(m['alignment_query_coverage'])<0.3 for m in models))
        self.assertTrue(all(m['evidence']=='tblastn_annotation_and_blastp' for m in models))
        self.assertEqual({r['status'] for r in table(out/'genome_query_summary.tsv') if r['query']=='q0002'},{'no_hit'})
        for protein in read_fasta(out/'sequences/q0001/protein.faa').values():
            self.assertEqual(len(protein),384)
            self.assertEqual(hashlib.sha256(protein.encode()).hexdigest(),'3c4596ae8dc9a468015fef0d36971a2a4b4a9271b89e546c79763b67f7400d4d')
        self.assertIn('local_alignment',(out/'report.html').read_text())

    def test_short_exon_seeds_are_not_filtered_by_full_query_coverage(self):
        from genecontext.sequence import CODE
        rng=random.Random(312)
        protein='M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=599))
        codons={aa:codon for codon,aa in CODE.items() if aa!='*'}
        coding=''.join(codons[aa] for aa in protein)+'TAA'
        sequence='N'*300;lines=['##gff-version 3'];parts=[]
        for i in range(12):
            exon=coding[i*150:(i+1)*150] if i<11 else coding[i*150:]
            start=len(sequence);sequence+=exon;parts.append((start,len(sequence)))
            if i<11:sequence+='GT'+''.join(rng.choices('ACGT',k=296))+'AG'
        sequence+='N'*300
        lines.extend([f'chr\ttest\tgene\t{parts[0][0]+1}\t{parts[-1][1]}\t.\t+\t.\tID=opaque',
                      f'chr\ttest\tmRNA\t{parts[0][0]+1}\t{parts[-1][1]}\t.\t+\t.\tID=opaque.t;Parent=opaque'])
        lines.extend(f'chr\ttest\tCDS\t{a+1}\t{b}\t.\t+\t0\tParent=opaque.t' for a,b in parts)
        fa=self.root/'short_exons.fna';write_fasta(fa,[('chr',sequence)])
        gff=self.root/'short_exons.gff3';gff.write_text('\n'.join(lines)+'\n')
        manifest=self.root/'short_exons.tsv';manifest.write_text('genome_id\tfasta\tgff3\nshort_exons\t'+str(fa)+'\t'+str(gff)+'\n')
        write_fasta(self.query,[('opaque_long_reference',protein)])
        out=self.run_search('short_exons',manifest)
        models=table(out/'models.tsv');self.assertTrue(any(r['source_model']=='opaque.t' for r in models))
        self.assertIn(protein,read_fasta(out/'sequences/q0001/protein.faa').values())
        # At least one seed is shorter than the candidate's 15%-coverage gate.
        hsps=(out/'genomes/short_exons/protein_discovery/tblastn.tsv').read_text().splitlines()
        self.assertTrue(any((int(x.split('\t')[5])-int(x.split('\t')[4])+1)/600<0.15 for x in hsps))

    def test_unannotated_nanog_refined_with_dense_seeds_and_native_strands(self):
        out=self.run_search('without_gff',self.unannotated)
        models=table(out/'models.tsv');self.assertTrue(models)
        self.assertEqual({m['genome_id'] for m in models},{'forward','reverse'})
        self.assertTrue(all(m['source']=='miniprot_rescue' and m['status']=='partial' for m in models))
        self.assertTrue(all(m['homology_support']=='local_alignment' and float(m['protein_bitscore'])>=80 for m in models))
        for gid in ['forward','reverse']:
            info=json.loads((out/'genomes'/gid/'protein_discovery/discovery.json').read_text())
            self.assertEqual(info['status'],'completed');self.assertGreater(info['coding_candidates'],0)
        self.assertTrue(all(r['source']=='tblastn' for r in table(out/'dna_regions.tsv')))
        self.assertIn('tblastn',json.loads((out/'manifest.json').read_text())['tools'])


if __name__=='__main__':unittest.main()
