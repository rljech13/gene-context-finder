"""Locus assembly, explicit review, and actual BLAST-backed integration."""
import json
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from unittest.mock import patch

from genecontext.cli import main, parser
from genecontext.comparison import table
from genecontext.locus_chaining import chain_hsps, group_loci, summarise
from genecontext.sequence import CODE, read_fasta, write_fasta, write_tsv
from genecontext.stage1 import apply_filters, approve, consolidate_refined, panel, refilter, search, verify


def h(qstart,qend,start,end,query='q1',strand='+',score=100,identity=.8,hid=None):
    return dict(hsp_id=hid or f'{query}_{qstart}_{start}',query=query,contig='chr',strand=strand,
        start=start,end=end,qstart=qstart,qend=qend,alignment_length=qend-qstart,qlen=100,
        slen=100000,bitscore=score,identity=identity,evalue=1e-20)


class Chains(unittest.TestCase):
    def test_spliced_forward_fragments_are_one_chain(self):
        chains=chain_hsps([h(0,50,1000,1150),h(50,100,1500,1650)])
        self.assertEqual(len(chains),1);self.assertEqual(chains[0]['coverage'],1)

    def test_reverse_chain_follows_transcription_not_coordinate_sort(self):
        chains=chain_hsps([h(0,50,1500,1650,strand='-'),h(50,100,1000,1150,strand='-')])
        self.assertEqual(len(chains),1);self.assertEqual(chains[0]['start'],1000)

    def test_adjacent_full_copies_are_never_proximity_merged(self):
        chains=chain_hsps([h(0,100,1000,1300),h(0,100,1400,1700)])
        self.assertEqual(len(chains),2);self.assertEqual(len(group_loci(chains)),2)

    def test_large_intron_requires_explicit_change(self):
        values=[h(0,50,1000,1150),h(50,100,3000,3150)]
        self.assertEqual(len(chain_hsps(values,max_intron=1000)),2)
        self.assertEqual(len(chain_hsps(values,max_intron=3000)),1)

    def test_query_resets_different_contigs_and_strands_do_not_chain(self):
        a,b=h(50,100,1000,1150),h(0,50,1500,1650)
        self.assertEqual(len(chain_hsps([a,b])),2)
        b['qstart'],b['qend'],b['contig']=50,100,'other'
        self.assertEqual(len(chain_hsps([h(0,50,1000,1150),b])),2)

    def test_overlapping_query_bases_are_not_double_counted(self):
        value=summarise([h(0,60,1000,1180),h(50,100,1150,1300)])
        self.assertEqual(value['coverage'],1);self.assertAlmostEqual(value['identity'],.8)

    def test_reference_isoforms_merge_on_shared_aligned_blocks(self):
        chains=chain_hsps([h(0,100,1000,1300),h(0,90,1015,1285,query='q2')])
        self.assertEqual(len(group_loci(chains)),1)

    def test_large_envelope_does_not_bridge_independent_copies(self):
        a=dict(contig='chr',strand='+',blocks=[(100,200)],start=100,end=200,bitscore=100,chain_id='a')
        b=dict(a,blocks=[(500,600)],start=500,end=600,chain_id='b')
        bridge=dict(a,blocks=[(100,200),(500,600)],start=100,end=600,bitscore=300,chain_id='bridge')
        clusters=group_loci([a,b,bridge]);self.assertEqual(len(clusters),2)

    def test_thresholds_must_pass_in_one_chain(self):
        queries={'a':dict(sequence='M'*100),'b':dict(sequence='M'*100)}
        locus=dict(genome_id='g',gene='TEST',contig='chr',start=100,locus_id='test',chains=[dict(chain_id='a',query='a',coverage=.9,identity=.2,evalue=1e-20,bitscore=100),
                          dict(chain_id='b',query='b',coverage=.2,identity=.9,evalue=1e-20,bitscore=100)],flags=[])
        apply_filters([locus],dict(min_query_coverage=.5,min_identity=.5,max_evalue=1e-5),queries)
        self.assertFalse(locus['pass_filters'])

    def test_model_anchors_split_a_chain_that_crosses_tandem_copies(self):
        values=[h(0,50,1000,1150),h(50,100,3000,3150)]
        chain=chain_hsps(values)[0]
        models=[dict(model_id='m1',query='q1',start=1000,end=1300,cds=[[1000,1300,0]],coverage=1,identity=1,rank=1),
                dict(model_id='m2',query='q1',start=2850,end=3150,cds=[[2850,3150,0]],coverage=1,identity=1,rank=2)]
        locus=dict(locus_id='old',genome_id='g',contig='chr',strand='+',chains=[chain],models=models,contig_length=100000)
        queries={'q1':dict(sequence='M'*100,members=[dict(query='variant',species='species')])}
        config=dict(gene='TEST',locus_shared_fraction=.5,window=50000)
        output=consolidate_refined([locus],queries,config)
        self.assertEqual(len(output),2)
        self.assertEqual([(l['start'],l['end']) for l in output],[(1000,1300),(2850,3150)])
        self.assertTrue(all('seed_chain_split_between_anchors' in l['flags'] for l in output))
        self.assertEqual({h['hsp_id'] for l in output for c in l['chains'] for h in c['hsps']}, {v['hsp_id'] for v in values})


@unittest.skipUnless(all(shutil.which(x) for x in ['samtools','tblastn','makeblastdb']),'requires actual BLAST and samtools')
class Stage1Integration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.root=Path(cls.tmp.name)
        rng=random.Random(316)
        amino='ACDEFGHIKLMNPQRSTVWY'
        protein='M'+''.join(rng.choice(amino) for _ in range(199))
        codons={aa:next(c for c,a in CODE.items() if a==aa) for aa in amino}
        coding=''.join(codons[a] for a in protein)
        cls.coding=coding+'TAA'
        background=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
        gene=coding[:300]+'GT'+background(996)+'AG'+coding[300:]+'TAA'
        # Two complete tandem copies plus a gene-free second genome and one
        # missing FASTA: all outcomes must survive in the review checkpoint.
        write_fasta(cls.root/'genome.fna',[('chr',background(1000)+gene+background(1000)+gene+background(1000))])
        write_fasta(cls.root/'empty.fna',[('chr',background(10000))])
        write_fasta(cls.root/'queries.faa',[('isoform_a',protein),('same_protein_other_transcript',protein),('isoform_b',protein[:-20])])
        write_tsv(cls.root/'queries.tsv',[
            dict(query='isoform_a',species='Species one'),
            dict(query='same_protein_other_transcript',species='Species one'),
            dict(query='isoform_b',species='Species two')],['query','species'])
        write_tsv(cls.root/'genomes.tsv',[
            dict(genome_id='copies',fasta='genome.fna'),dict(genome_id='empty',fasta='empty.fna'),
            dict(genome_id='broken',fasta='missing.fna')],['genome_id','fasta'])
        cls.output=cls.root/'result'
        cls.argv=['stage1','--gene','TEST','--query',str(cls.root/'queries.faa'),'--query-type','protein',
            '--query-metadata',str(cls.root/'queries.tsv'),'--genomes',str(cls.root/'genomes.tsv'),
            '--out',str(cls.output),'--refine','none','--threads','2']
        cls.code=search(parser().parse_args(cls.argv))

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def edit_decisions(self,name='decisions.tsv'):
        rows=table(self.output/'loci.decisions.tsv')
        for row in rows:
            row['decision']='inspect' if row['locus_id'] else 'unresolved'
            row['notes']='Integration test, not a biological approval'
        path=self.root/name;write_tsv(path,rows,list(rows[0]));return path,rows

    def test_real_blast_recovers_two_spliced_copies_and_collapses_isoforms(self):
        _,data=verify(self.output)
        own=[l for l in data['loci'] if l['genome_id']=='copies' and l['pass_filters']]
        self.assertEqual(len(own),2)
        for locus in own:
            self.assertEqual(locus['query_variants'],3)
            self.assertEqual(locus['unique_sequences'],2)
            self.assertEqual(locus['support_species'],2)
            self.assertGreater(locus['coverage'],.95)

    def test_nohit_and_failed_genomes_are_explicit_and_no_automatic_decisions(self):
        self.assertEqual(self.code,1)
        summary={r['genome_id']:r for r in table(self.output/'genomes.summary.tsv')}
        self.assertEqual(summary['empty']['status'],'no_candidate')
        self.assertEqual(summary['broken']['status'],'error')
        self.assertTrue(all(not r['decision'] for r in table(self.output/'loci.decisions.tsv')))

    def test_refilter_uses_saved_evidence_without_subprocesses(self):
        destination=self.root/'filtered'
        with patch('subprocess.run',side_effect=AssertionError('No tools should run')):
            refilter(parser().parse_args(['stage1-filter','--results',str(self.output),'--out',str(destination),'--min-identity','1']))
        self.assertEqual(len(table(destination/'loci.tsv')),len(table(self.output/'loci.tsv')))
        self.assertNotEqual(verify(destination)[1]['review_id'],verify(self.output)[1]['review_id'])

    def test_pending_decisions_block_before_output_creation(self):
        out=self.root/'pending_out'
        args=parser().parse_args(['stage1-select','--results',str(self.output),'--decisions',str(self.output/'loci.decisions.tsv'),'--out',str(out)])
        with self.assertRaisesRegex(ValueError,'(?i)pending'):approve(args)
        self.assertFalse(out.exists())

    def test_manual_windows_and_split_anchors_are_exported(self):
        path,rows=self.edit_decisions()
        row=next(r for r in rows if r['locus_id']);start,end=map(int,row['anchor_intervals'].split(':'))
        row['anchor_intervals']=f'{start}:{start+100};{end-100}:{end}'
        write_tsv(path,rows,list(rows[0]));out=self.root/'approved'
        approve(parser().parse_args(['stage1-select','--results',str(self.output),'--decisions',str(path),'--out',str(out)]))
        selected=table(out/'approved_loci.tsv');self.assertEqual(len(selected),3)
        self.assertTrue(any(r['region_id'].endswith('__part2') for r in selected))
        self.assertTrue((out/'synteny_windows.bed').exists())

    def test_unknown_or_cross_checkpoint_decisions_are_rejected(self):
        path,rows=self.edit_decisions('wrong.tsv');rows[0]['review_id']='wrong';write_tsv(path,rows,list(rows[0]))
        with self.assertRaisesRegex(ValueError,'checkpoint'):
            approve(parser().parse_args(['stage1-select','--results',str(self.output),'--decisions',str(path),'--out',str(self.root/'wrongout')]))

    def test_editing_decisions_is_allowed_but_changing_evidence_is_detected(self):
        copy=self.root/'tamper_check';shutil.copytree(self.output,copy)
        (copy/'loci.decisions.tsv').write_text('An editable working decision file\n')
        verify(copy)
        with open(copy/'loci.tsv','a') as handle:handle.write('changed evidence\n')
        with self.assertRaisesRegex(ValueError,'evidence changed'):verify(copy)

    def test_cds_and_dna_queries_both_search_the_supplied_genome(self):
        query=self.root/'coding.fna';write_fasta(query,[('coding_reference',self.coding)])
        manifest=self.root/'one.tsv';write_tsv(manifest,[dict(genome_id='copies',fasta='genome.fna')],['genome_id','fasta'])
        for kind in ['cds','dna']:
            output=self.root/('route_'+kind)
            args=parser().parse_args(['stage1','--gene','TEST','--query',str(query),'--query-type',kind,
                '--genomes',str(manifest),'--out',str(output),'--refine','none','--threads','2'])
            self.assertEqual(search(args),0)
            _,data=verify(output)
            self.assertEqual(len([l for l in data['loci'] if l['pass_filters']]),2)
            self.assertTrue(any(c['argv'][0]==('tblastn' if kind=='cds' else 'blastn') for c in data['commands']))

    def test_parallel_genomes_match_serial_candidates_and_keep_failures(self):
        from genecontext.stage1_parallel import search as parallel
        args=parser().parse_args(self.argv);args.out=str(self.root/'parallel');args.jobs=2
        self.assertEqual(parallel(args),1)  # The deliberately missing genome remains an error.
        _,serial=verify(self.output);_,batched=verify(args.out)
        fields=['genome_id','locus_id','start','end','strand','coverage','identity','flags']
        signature=lambda data:sorted(tuple(str(l[k]) for k in fields) for l in data['loci'])
        self.assertEqual(signature(serial),signature(batched))
        self.assertEqual({g['genome_id']:g['status'] for g in serial['genomes']},
                         {g['genome_id']:g['status'] for g in batched['genomes']})
        self.assertTrue(all(not r['decision'] for r in table(Path(args.out)/'loci.decisions.tsv')))

    def test_merge_rejects_duplicate_genomes(self):
        from genecontext.stage1_parallel import merge_parts
        with self.assertRaisesRegex(ValueError,'Duplicate genome'):
            merge_parts([self.output,self.output],self.root/'invalid_merge',['copies','empty','broken'])


class QueryBundle(unittest.TestCase):
    def test_shipped_annotation_cds_translations_and_provenance(self):
        root=Path(__file__).resolve().parents[1]/'examples/queries/pou5_test'
        from genecontext.sequence import translate
        proteins,cds=read_fasta(root/'proteins.faa'),read_fasta(root/'cds.fna')
        rows=table(root/'queries.tsv')
        self.assertEqual(len(rows),16);self.assertEqual(set(proteins),set(cds))
        self.assertEqual(len({r['species'] for r in rows}),4)
        for row in rows:
            translated=translate(cds[row['query']]).removesuffix('*')
            if row['initiation_exception']=='True':self.assertEqual(translated[1:],proteins[row['query']][1:])
            else:self.assertEqual(translated,proteins[row['query']])
        for source in json.loads((root/'provenance.json').read_text())['sources']:
            from genecontext.sequence import sha256
            self.assertEqual(sha256(root/source['file']),source['sha256'])


@unittest.skipUnless(all(shutil.which(x) for x in ['samtools','tblastn','makeblastdb','miniprot']),'requires BLAST, samtools, miniprot')
class POU5ReferenceControls(unittest.TestCase):
    def test_annotated_reference_windows_recover_genes_and_keep_frog_copies_apart(self):
        panel_root=Path(__file__).resolve().parents[1]/'examples/queries/pou5_test'
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'out'
            args=parser().parse_args(['stage1','--gene','POU5','--query',str(panel_root/'proteins.faa'),
                '--query-type','protein','--query-metadata',str(panel_root/'queries.tsv'),
                '--genomes',str(panel_root/'controls.genomes.tsv'),'--out',str(root),'--threads','2','--refine','miniprot'])
            self.assertEqual(search(args),0)
            _,data=verify(root)
            self.assertTrue(all(g['control_overlap'] for g in data['genomes']))
            expected={'Homo_sapiens__POU5F1':1,'Danio_rerio__pou5f3':1,'Gallus_gallus__Pou5f3':1,
                      'Xenopus_tropicalis__pou5f3':2,'Xenopus_tropicalis__pou5f3.2':3,'Xenopus_tropicalis__pou5f3.3':2}
            for genome,count in expected.items():
                own=[l for l in data['loci'] if l['genome_id']==genome]
                self.assertEqual(len(own),count,genome)
                self.assertTrue(all(l['models'] for l in own))
            # Every retained BLAST seed occurs exactly once after consolidation.
            for genome in data['genomes']:
                raw=json.loads((root/'raw'/(genome['genome_id']+'.hsps.json')).read_text())
                dropped={h['hsp_id'] for h in raw['discarded_seeds']}
                expected_hsps={h['hsp_id'] for h in raw['hsps']}-dropped
                assigned=[h['hsp_id'] for l in data['loci'] if l['genome_id']==genome['genome_id'] for c in l['chains'] for h in c['hsps']]
                self.assertEqual(set(assigned),expected_hsps)
                self.assertEqual(len(assigned),len(set(assigned)))


if __name__=='__main__':unittest.main()
