"""Comparative controls: gene order, missing evidence, ambiguity and base accounting."""
import csv
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from genecontext.comparison import make_groups, orient_interval, select_candidates
from genecontext.context_compare import synteny_pairs
from genecontext.dna_compare import parse_lastz, summarize_blocks
from genecontext.sequence import CODE, read_fasta, read_gff, reverse_complement, translate, write_fasta


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class ComparisonUnits(unittest.TestCase):
    def test_group_mapping_requires_every_query_once(self):
        queries = [dict(id='q0001', name='one'), dict(id='q0002', name='two')]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'groups.tsv'
            path.write_text('query\tgroup\none\tgene\nq0002\tgene\n')
            self.assertEqual(make_groups(queries, path)[0]['queries'], ['q0001', 'q0002'])
            path.write_text('query\tgroup\none\tgene\n')
            with self.assertRaises(ValueError):
                make_groups(queries, path)

    def test_isoforms_collapse_but_separate_loci_are_retained(self):
        base = dict(query='q0001', genome_id='g', contig='c', strand='+', status='complete_predicted',
                    source='annotation', source_gene='one', protein_aa=100, coverage=1, identity=1)
        a = dict(base, model_id='a', start=100, end=403, cds_exons='[[100,403,0]]')
        b = dict(base, model_id='b', start=100, end=400, cds_exons='[[100,400,0]]')
        c = dict(base, model_id='c', source_gene='two', start=1000, end=1303, cds_exons='[[1000,1303,0]]')
        selected, candidates, summary = select_candidates([dict(id='q0001', name='gene', queries=['q0001'])],
                                                         [dict(genome_id='g')], [a, b, c], [])
        self.assertEqual(summary[0]['candidate_loci'], 2)
        self.assertEqual(summary[0]['status'], 'ranked_multiple_loci')
        self.assertEqual(len(candidates), 3)
        self.assertEqual(selected[0]['model_id'], 'a')
        # A strong incomplete target must not lose to a weak complete paralog.
        strong_partial=dict(c,status='partial',coverage=0.9,identity=0.95)
        weak_complete=dict(a,coverage=0.6,identity=0.4)
        chosen,_,_=select_candidates([dict(id='q0001',name='gene',queries=['q0001'])],
                                    [dict(genome_id='g')],[weak_complete,strong_partial],[])
        self.assertEqual(chosen[0]['model_id'],'c')

    def test_synteny_uses_order_and_unique_anchors(self):
        loci = [dict(locus_id=x, group='q', context_status='provided_gff3') for x in ['a','b']]
        features = [dict(locus_id=lid, homology_group=h, relative_start=i*100, relative_strand='+')
                    for lid, groups in [('a',['H1','H2','H3']),('b',['H3','H2','H1'])]
                    for i,h in enumerate(groups)]
        result = synteny_pairs(loci, features)[0]
        self.assertEqual(result['unique_anchors'], 3)
        self.assertEqual(result['order_concordance'], 0)
        self.assertEqual(result['status'], 'different_local_order_or_orientation')
        loci[1]['context_status'] = 'no_annotation_evidence'
        self.assertEqual(synteny_pairs(loci, features)[0]['status'], 'context_evidence_missing')
        features.extend([dict(features[0], relative_start=900), dict(features[1], relative_start=1000)])
        self.assertLess(synteny_pairs(loci, features)[0]['unique_anchors'], 3)

    def test_dna_gaps_ambiguity_and_overlaps_do_not_inflate_coverage(self):
        a,b = 'ACGTNACGT','ACGTTACGT'
        block = dict(start_a=0,end_a=9,start_b=0,end_b=9,first_a=0,first_b=0,
                     step_a=1,step_b=1,text_a=a,text_b=b,score=100)
        stats, blocks, columns = summarize_blocks([block,dict(block,score=50)],a,b,4,0.5)
        self.assertEqual(stats['paired_bases'], 8)
        self.assertEqual(stats['identity'], 1)
        self.assertEqual(stats['coverage_a'], 1)
        self.assertEqual(stats['coverage_b'], 8/9)
        self.assertEqual(blocks[1]['status'],'overlap_with_higher_scoring_block')
        self.assertEqual(len(columns),8)

    def test_lastz_reverse_coordinates_roundtrip(self):
        seq='ACGTTGCAACGA'
        other=reverse_complement(seq)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'lastz.tsv'
            path.write_text('#score\n100\t0\t12\t0\t12\t-\t'+seq+'\t'+seq+'\n')
            stats,blocks,_=summarize_blocks(parse_lastz(path),seq,other,4,0.5)
            self.assertEqual(stats['identity'],1)
            self.assertEqual(blocks[0]['relative_strand'],'-')

    def test_minus_window_and_lowercase_mask_preservation(self):
        row=dict(strand='-',window_start=100,window_end=900)
        self.assertEqual(orient_interval(row,200,300),(600,700))
        self.assertEqual(reverse_complement('AaCcGgTtNn'),'nNaAcCgGtT')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mask.fa';path.write_text('>x\nAaCc\n')
            self.assertEqual(read_fasta(path)['x'],'AACC')
            self.assertEqual(read_fasta(path,preserve_case=True)['x'],'AaCc')


REQUIRED=['samtools','miniprot','blastn','blastp','blastx','makeblastdb','mafft']


@unittest.skipUnless(all(shutil.which(p) for p in REQUIRED), 'requires actual sequence tools')
class ComparativeIntegration(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        rng=random.Random(812)
        codons={}
        for codon,aa in CODE.items():codons.setdefault(aa,[]).append(codon)
        self.proteins={key:'M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=179)) for key in ['A','TARGET','B','absent']}
        self.coding={key:''.join(rng.choice(codons[aa]) for aa in p)+'TAA' for key,p in self.proteins.items()}
        self.manifest=self.root/'genomes.tsv'
        rows=[]
        for gid,order,reverse,annotated in [('first',['A','TARGET','B'],False,True),
                                            ('rearranged',['B','TARGET','A'],True,True),
                                            ('unannotated',['A','TARGET','B'],True,False)]:
            seq=''.join(rng.choices('ACGT',k=800)).lower();features=[]
            for gene in order:
                start=len(seq);seq+=self.coding[gene];end=len(seq)
                features.append([gene,start,end,'+'])
                seq+=''.join(rng.choices('ACGT',k=600))
            if reverse:
                n=len(seq);seq=reverse_complement(seq)
                features=[[g,n-e,n-s,'-'] for g,s,e,_ in features]
            path=self.root/(gid+'.fna.gz')
            with gzip.open(path,'wt') as h:h.write('>chr_test\n'+seq+'\n')
            gff=self.root/(gid+'.gff3');lines=['##gff-version 3']
            for gene,start,end,strand in features:
                # Different annotation names between genomes must not prevent sequence-based matching.
                label=gid+'_'+gene
                lines += [f'chr_test\ttest\tgene\t{start+1}\t{end}\t.\t{strand}\t.\tID={label};Name={label}',
                          f'chr_test\ttest\tmRNA\t{start+1}\t{end}\t.\t{strand}\t.\tID={label}.t;Parent={label}',
                          f'chr_test\ttest\tCDS\t{start+1}\t{end}\t.\t{strand}\t0\tParent={label}.t']
            gff.write_text('\n'.join(lines)+'\n')
            rows.append((gid,str(path),str(gff) if annotated else ''))
        self.rows=rows
        self.manifest.write_text('genome_id\tfasta\tgff3\n'+'\n'.join('\t'.join(r) for r in rows)+'\n')
        self.query=self.root/'query.faa'
        write_fasta(self.query,[('target',self.proteins['TARGET']),('negative',self.proteins['absent'])])

    def tearDown(self):self.temp.cleanup()

    def execute(self,command,expected=0):
        # Preserve v0.2 annotation-first regression controls explicitly.
        if command[0] in {'analyze', 'compare'} and '--context-mode' not in command:
            command = [*command, '--context-mode', 'annotation', '--context-reference', 'none']
        run=subprocess.run([sys.executable,'-m','genecontext',*map(str,command)],capture_output=True,text=True)
        self.assertEqual(run.returncode,expected,run.stdout+run.stderr)
        return run

    def test_full_analysis_synteny_dna_proteins_and_unannotated_context(self):
        output=self.root/'analysis'
        self.execute(['analyze','--query',self.query,'--query-type','protein','--genomes',self.manifest,
                      '--out',output,'--threads','1','--window','5000'])
        comparison=output/'comparison'
        manifest=json.loads((comparison/'manifest.json').read_text())
        self.assertEqual(manifest['status'],'completed',manifest.get('errors'))
        self.assertEqual(len(table(comparison/'summary.tsv')),6)
        self.assertEqual(len(table(comparison/'loci.tsv')),3)
        neighbors=table(comparison/'neighbors.tsv')
        inferred=[r for r in neighbors if r['genome_id']=='unannotated']
        self.assertEqual(len(inferred),2,inferred)
        self.assertTrue(all(r['evidence']=='protein_supported_prediction' for r in inferred))
        synteny=table(comparison/'synteny.tsv')
        loci={r['locus_id']:r for r in table(comparison/'loci.tsv')}
        ordered={(loci[r['locus_a']]['genome_id'],loci[r['locus_b']]['genome_id']):r for r in synteny}
        self.assertEqual(ordered['first','rearranged']['status'],'different_local_order_or_orientation')
        self.assertEqual(ordered['first','unannotated']['status'],'conserved_local_order')
        self.assertEqual(len(table(comparison/'dna_pairs.tsv')),3)
        self.assertTrue(all(r['status']=='aligned_blocks' for r in table(comparison/'dna_pairs.tsv')))
        self.assertTrue(all(0 < float(r['coverage_a']) <= 1 and 0 < float(r['coverage_b']) <= 1 for r in table(comparison/'dna_pairs.tsv')))
        self.assertTrue(all(r['left_clipped']=='True' for r in loci.values()))
        self.assertTrue(any(int(r['softmasked_bases'])>0 for r in loci.values()))
        self.assertEqual(next(r['target_name'] for r in loci.values() if r['genome_id']=='first'),'first_TARGET')
        self.assertEqual({r['category'] for r in table(comparison/'dna_categories.tsv')},{'coding','upstream','downstream'})
        aligned=read_fasta(comparison/'queries/q0001/proteins.aligned.faa')
        self.assertEqual(len(aligned),4)
        codons=read_fasta(comparison/'queries/q0001/cds.codon.aligned.fna')
        for key,seq in codons.items():self.assertEqual(translate(seq.replace('-','')),self.proteins['TARGET'])
        models,_=read_gff(comparison/'context_predictions/unannotated.gff3')
        self.assertEqual(len(models),2)
        self.assertTrue(all(m.strand=='-' for m in models))
        native=read_fasta(self.rows[2][1])['chr_test']
        for model in models:
            dna=''.join(reverse_complement(native[a:b]) for a,b,_ in sorted(model.cds,reverse=True))
            self.assertIn(translate(dna).removesuffix('*'),{self.proteins['A'],self.proteins['B']})
        for path in comparison.rglob('*.svg'):ET.parse(path)
        self.assertTrue((output/'report.html').exists())

    def test_no_annotation_evidence_and_grouped_queries_are_explicit(self):
        manifest=self.root/'unannotated.tsv'
        manifest.write_text('genome_id\tfasta\tgff3\n'+self.rows[2][0]+'\t'+self.rows[2][1]+'\n')
        write_fasta(self.query,[('ref_one',self.proteins['TARGET']),('ref_two',self.proteins['TARGET'])])
        groups=self.root/'groups.tsv';groups.write_text('query\tgroup\nref_one\tMY_GENE\nref_two\tMY_GENE\n')
        output=self.root/'grouped'
        self.execute(['analyze','--query',self.query,'--query-type','protein','--genomes',manifest,
                      '--out',output,'--query-groups',groups,'--threads','1'])
        loci=table(output/'comparison/loci.tsv')
        self.assertEqual(len(loci),1)
        self.assertEqual(loci[0]['context_status'],'no_annotation_evidence')
        self.assertEqual(table(output/'comparison/summary.tsv')[0]['candidate_loci'],'1')
        self.assertEqual(len(table(output/'comparison/candidates.tsv')),2)
        self.assertIn('MY_GENE',(output/'comparison/report.html').read_text())

    def test_noncoding_query_remains_dna_and_compare_preserves_source(self):
        query=self.root/'region.fna'
        native=read_fasta(self.rows[0][1])['chr_test']
        write_fasta(query,[('masked_region',native[:500]),('mixed_region',native[600:1100])])
        manifest=self.root/'dna.tsv'
        manifest.write_text('genome_id\tfasta\tgff3\nfirst\t'+self.rows[0][1]+'\nsecond\t'+self.rows[0][1]+'\n')
        output=self.root/'dna_analysis'
        self.execute(['analyze','--query',query,'--query-type','dna','--genomes',manifest,
                      '--out',output,'--threads','1','--window','300'])
        comparison=output/'comparison'
        self.assertEqual({r['status'] for r in table(comparison/'loci.tsv')},{'dna_region'})
        pairs={r['group']:r for r in table(comparison/'dna_pairs.tsv')}
        self.assertEqual(pairs['q0001']['status'],'no_passing_blocks')
        self.assertEqual(pairs['q0001']['identity'],'')
        self.assertEqual(pairs['q0002']['status'],'aligned_blocks')
        original=(output/'search/manifest.json').read_bytes()
        self.execute(['compare','--results',output/'search','--out',self.root/'again','--threads','1','--window','200'])
        self.assertEqual((output/'search/manifest.json').read_bytes(),original)

    def test_external_context_references_and_failed_genome_remain_visible(self):
        manifest=self.root/'broken_and_raw.tsv'
        manifest.write_text('genome_id\tfasta\tgff3\nbroken\tmissing.fna\nraw\t'+self.rows[2][1]+'\n')
        write_fasta(self.query,[('arbitrary_target',self.proteins['TARGET'])])
        refs=self.root/'references.faa'
        write_fasta(refs,[('reference_alpha',self.proteins['A']),('reference_beta',self.proteins['B'])])
        output=self.root/'external'
        self.execute(['analyze','--query',self.query,'--query-type','protein','--genomes',manifest,
                      '--out',output,'--context-proteins',refs,'--threads','1'],expected=1)
        summary=table(output/'comparison/summary.tsv')
        self.assertEqual(summary[0]['status'],'source_search_error')
        self.assertEqual(summary[1]['status'],'single_candidate_locus')
        self.assertEqual(len(table(output/'comparison/neighbors.tsv')),2)
        self.assertEqual(json.loads((output/'comparison/manifest.json').read_text())['status'],'completed_with_errors')
        self.assertIn('arbitrary_target',(output/'report.html').read_text())

    @unittest.skipUnless(shutil.which('lastz'),'optional LASTZ executable required')
    def test_actual_lastz_reverse_alignment(self):
        from genecontext.dna_compare import parse_lastz
        from genecontext.tools import Runner
        rng=random.Random(19)
        sequence=''.join(rng.choices('ACGT',k=900))
        a,b=self.root/'a.fa',self.root/'b.fa'
        write_fasta(a,[('one',sequence)])
        write_fasta(b,[('two',reverse_complement(sequence))])
        output=self.root/'lastz.tsv'
        Runner(self.root).run(['lastz',a,b,'--strand=both','--format=general:score,zstart1,end1,zstart2+,end2+,strand2,text1,text2'],output)
        stats,blocks,_=summarize_blocks(parse_lastz(output),sequence,reverse_complement(sequence))
        self.assertEqual(stats['coverage_a'],1)
        self.assertEqual(stats['coverage_b'],1)
        self.assertEqual(stats['identity'],1)
        self.assertEqual({r['relative_strand'] for r in blocks if r['status']=='accepted'},{'-'})

    def test_clean_partial_protein_option_and_codon_roundtrip(self):
        from genecontext.cli import parser
        from genecontext.comparison_report import protein_comparison
        from genecontext.tools import Runner
        source=self.root/'saved_search';seqdir=source/'sequences/q0001';seqdir.mkdir(parents=True)
        fragment=self.coding['TARGET'][30:-3]
        write_fasta(seqdir/'cds.fna',[('partial_model',fragment)])
        write_fasta(seqdir/'protein.faa',[('partial_model',translate(fragment))])
        output=self.root/'partial_output';output.mkdir()
        args=parser().parse_args(['compare','--results',str(source),'--out',str(output),'--threads','1','--include-partial'])
        groups=[dict(id='q0001',name='target',queries=['q0001'])]
        queries=[dict(id='q0001',protein=self.proteins['TARGET'],coding='')]
        selected=[dict(group='q0001',query='q0001',model_id='partial_model',locus_id='L000001',status='partial')]
        errors=[]
        summary=protein_comparison(args,source,output,groups,queries,selected,Runner(output),errors)
        self.assertFalse(errors)
        self.assertEqual(summary[0]['candidate_proteins'],1)
        codons=read_fasta(output/'queries/q0001/cds.codon.aligned.fna')
        self.assertEqual(codons['L000001'].replace('-',''),fragment)
        self.assertEqual(translate(fragment),self.proteins['TARGET'][10:])


if __name__=='__main__':unittest.main()
