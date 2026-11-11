"""STOP 2: approval provenance, coordinates, real BLAST and prediction adapter."""
import copy
import gzip
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.sequence import Model,sha256,write_fasta,write_tsv
from genecontext.stage2 import selected_regions,extract_intervals,model_record,assign_groups,run
from genecontext.stage2_db import prepare,verify_db


class DiagramSourceEvidence(unittest.TestCase):
    def fixture(self):
        models = [dict(model_id=source+genome, source=source, genome_id=genome,
                       gene_group=genome, overlaps_anchor=False)
                  for source in ['annotation', 'tiberius'] for genome in ['a', 'b']]
        return dict(models=models, neighbor_pairs=[dict(model_a='annotationa',
                    model_b='annotationb', gene_a='a', gene_b='b', bitscore=100)],
                    provenance=dict(parameters=dict(identification_coverage=.5,
                    identification_identity=.3, competing_ratio=.95)))

    def test_source_links_use_saved_prediction_proteins_not_annotation_winner(self):
        from genecontext.stage2_report import source_pairs
        data=self.fixture();before=copy.deepcopy(data)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'neighbor_support').mkdir()
            write_fasta(root/'proteins.faa',[(m['model_id'],'M'*10) for m in data['models']])
            (root/'neighbor_support/blast.tsv').write_text(''.join(
                f'{s}a\t{s}b\t80\t10\t1\t10\t1\t10\t1e-10\t{score}\t10\t10\n'
                for s,score in [('annotation',100),('tiberius',80)]))
            pairs=source_pairs(root,data)
            self.assertEqual(pairs['tiberius'][0]['model_a'],'tiberiusa')
            self.assertEqual(pairs['tiberius'][0]['bitscore'],80)
            self.assertEqual(pairs['annotation'][0]['bitscore'],100)
            self.assertEqual(data,before)

    def test_missing_raw_hits_never_borrow_annotation_support_for_tiberius(self):
        from genecontext.stage2_report import source_pairs
        with tempfile.TemporaryDirectory() as d:
            pairs=source_pairs(Path(d),self.fixture())
            self.assertEqual(pairs['tiberius'],[])
            self.assertEqual(len(pairs['annotation']),1)


class SelectionAndCoordinates(unittest.TestCase):
    def data(self):
        return dict(review_id='one',config=dict(gene='TEST',window=50),
            genomes=[dict(genome_id='a',status='success'),dict(genome_id='b',status='success')],
            loci=[dict(genome_id=g,locus_id=g,gene='TEST',contig='chr',strand='+',start=100,end=200,
                       contig_length=1000,pass_filters=True,boundary_source='annotation') for g in ['a','b']])
    def test_unsaved_checks_are_never_approved_and_stale_ids_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'draft.json';value=dict(schema='stage1_pair_draft_v1',review_id='one',pairs=[
                dict(genome_id=g,selected=[g],committed=g=='a',action='inspect',reason='sequence_support',notes='') for g in ['a','b']])
            path.write_text(json.dumps(value));args=SimpleNamespace(selection=None,decisions=None,draft=str(path),committed_only=False)
            with self.assertRaisesRegex(ValueError,'unsaved'):selected_regions(self.data(),args)
            args.committed_only=True;regions,receipt,_=selected_regions(self.data(),args)
            self.assertEqual([r['locus_id'] for r in regions],['a']);self.assertEqual(receipt['omitted_unsaved_genomes'],['b'])
            value['review_id']='wrong';path.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError,'Wrong'):selected_regions(self.data(),args)
    def test_streaming_gzip_overlaps_boundary_contigs_and_duplicate_names(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'g.fa.gz'
            with gzip.open(path,'wt') as out:out.write('>chr\nACGTac\ngtACGT\n>other\nNNNN\n')
            seq,length=extract_intervals(path,dict(a=('chr',2,10),b=('chr',4,8),c=('other',1,4)))
            self.assertEqual(seq,dict(a='GTacgtAC',b='acgt',c='NNN'));self.assertEqual(length['chr'],12)
            with self.assertRaisesRegex(ValueError,'outside'):extract_intervals(path,dict(bad=('other',0,5)))
            with gzip.open(path,'wt') as out:out.write('>chr\nAAA\n>chr\nAAA\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):extract_intervals(path,{})
    def test_native_reverse_cds_coordinates_and_fused_predictions(self):
        region=dict(region_id='r',locus_id='l',genome_id='g',contig='chr',window_start=1000,window_end=2000,anchor_start=1100,anchor_end=1400,strand='-')
        m=Model('p','r','-',cds=[(300,400,0),(100,200,2)],stops=[(97,100)],gene='p')
        rec=model_record(m,region,'tiberius','ATG','M','partial',[])
        self.assertEqual(rec['start'],1097);self.assertEqual(rec['end'],1400);self.assertEqual(rec['cds_blocks'],[[1300,1400,0],[1100,1200,2]])
        a=dict(rec,model_id='a',source='annotation',gene_group='A',cds_blocks=[[1100,1200,0]])
        b=dict(a,model_id='b',gene_group='B',cds_blocks=[[1300,1400,0]])
        assign_groups([a,b,rec]);self.assertEqual(set(rec['overlapping_annotation_groups']),{'A','B'});self.assertIn('overlaps_multiple_annotated_genes',rec['issues'])
        self.assertNotIn(rec['gene_group'],{'A','B'})


@unittest.skipUnless(shutil.which('blastp') and shutil.which('makeblastdb'),'requires BLAST+')
class RealBlastStage2(unittest.TestCase):
    def test_real_blast_stop2_with_injected_predictor_and_no_annotation(self):
        # Only the expensive neural predictor is injected. Extraction, CDS,
        # BLASTP database, matching, checkpoint and HTML are real.
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);query='M'+'ACDEFGHIKLMNPQRSTVWY'*8
            codons=dict(A='GCT',C='TGT',D='GAT',E='GAA',F='TTT',G='GGT',H='CAT',I='ATT',K='AAA',L='CTG',M='ATG',N='AAT',P='CCT',Q='CAA',R='CGT',S='TCT',T='ACT',V='GTT',W='TGG',Y='TAT')
            cds=''.join(codons[a] for a in query)+'TAA';dna='C'*1000+cds+'C'*1000
            write_fasta(root/'g.fna',[('chr',dna)]);write_fasta(root/'refs.faa',[('sp|A|ONE protein OS=species OX=9606 GN=TEST',query),('sp|B|TWO protein OS=species OX=9606 GN=OTHER',query)])
            prepare(SimpleNamespace(out=str(root/'db'),fasta=str(root/'refs.faa')));verify_db(root/'db')
            source=root/'source';source.mkdir()
            region_locus=dict(genome_id='g',locus_id='l',gene='TEST',contig='chr',strand='+',start=1000,end=1000+len(cds),contig_length=len(dna),pass_filters=True,boundary_source='test')
            data=dict(review_id='one',config=dict(gene='TEST',window=500,query_type='protein'),queries=[dict(id='q',sequence=query,members=[dict(query='reference')])],genomes=[dict(genome_id='g',status='success',fasta_sha256=sha256(root/'g.fna'))],loci=[region_locus])
            (source/'evidence.json').write_text(json.dumps(data));write_tsv(root/'genomes.tsv',[dict(genome_id='g',fasta=str(root/'g.fna'))],['genome_id','fasta'])
            write_tsv(root/'decisions.tsv',[dict(review_id='one',genome_id='g',locus_id='l',decision='inspect',reason='sequence_support',notes='unit control')],['review_id','genome_id','locus_id','decision','reason','notes'])
            args=parser().parse_args(['stage2','--results',str(source),'--decisions',str(root/'decisions.tsv'),'--genomes',str(root/'genomes.tsv'),'--protein-db',str(root/'db'),'--out',str(root/'result')])
            def predictor(args,out,rows,windows,runner,errors):
                rid=rows[0]['locus_id'];return {rid:[Model('p',rid,'+',cds=[(500,500+len(cds),0)],gene='p')]},dict(test_double=True)
            with patch('genecontext.stage2.verify',return_value=(source,data)),patch('genecontext.stage2.validate_prediction_inputs'),patch('genecontext.tiberius_backend.predict_batches',side_effect=predictor):
                self.assertEqual(run(args),0)
            result=json.loads((root/'result/evidence.json').read_text());self.assertEqual(len(result['models']),1)
            m=result['models'][0];self.assertEqual(m['coding_status'],'complete_predicted');self.assertEqual(m['identification']['status'],'ambiguous_family_assignment');self.assertEqual(m['identification']['proposed_name'],'');self.assertEqual(len(m['identification']['near_best']),2)
            self.assertTrue(m['target_support']);self.assertTrue((root/'result/report.html').exists());self.assertFalse((root/'result/protein_msa.faa').exists())
            seal=json.loads((root/'result/seal.json').read_text());self.assertNotIn('models.decisions.tsv',seal)
            self.assertTrue(all(sha256(root/'result'/name)==value for name,value in seal.items()))
            # Repeat only identification; none of the prediction or sequence
            # outputs may be recomputed, changed, or auto-approved.
            context_feature=dict(feature_id='rna',region_id='l',genome_id='g',original_gene='rna',name='TEST_RNA',kind='gene',biotype='lncRNA',contig='chr',start=600,end=650,strand='+')
            (root/'result/annotation_context.json').write_text(json.dumps(dict(schema='stage2_annotation_context_v1',source_evidence_sha256=sha256(root/'result/evidence.json'),features=[context_feature],sources=[])))
            seal['annotation_context.json']=sha256(root/'result/annotation_context.json');(root/'result/seal.json').write_text(json.dumps(seal))
            # Even an empty original-annotation supplement has an evidence
            # binding that must be refreshed when identification is replayed.
            (root/'result/annotation_quality.json').write_text(json.dumps(dict(schema='stage2_annotation_quality_v1',source_evidence_sha256=sha256(root/'result/evidence.json'),sources=[],models={})))
            seal['annotation_quality.json']=sha256(root/'result/annotation_quality.json');(root/'result/seal.json').write_text(json.dumps(seal))
            from genecontext.stage2 import reidentify
            replay=parser().parse_args(['stage2-reidentify','--results',str(root/'result'),'--protein-db',str(root/'db'),'--out',str(root/'reidentified')])
            self.assertEqual(reidentify(replay),0)
            updated=json.loads((root/'reidentified/evidence.json').read_text())
            self.assertNotEqual(updated['review_id'],result['review_id'])
            self.assertEqual(updated['annotation_features'],[context_feature])
            quality=json.loads((root/'reidentified/annotation_quality.json').read_text())
            self.assertEqual(quality['source_evidence_sha256'],sha256(root/'reidentified/evidence.json'))
            self.assertFalse((root/'reidentified/annotation_context.json').exists(),'The old binding must not survive a new evidence identity')
            for name in ['windows.fna','cds.fna','proteins.faa','models.gff3']:
                self.assertEqual(sha256(root/'result'/name),sha256(root/'reidentified'/name))
            self.assertEqual(updated['models'][0]['identification']['status'],'ambiguous_family_assignment')

    def test_symbol_case_is_not_a_competing_gene_and_accessions_stay_distinct(self):
        from genecontext.stage2 import identify
        from genecontext.tools import Runner
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);seq='MACDEFGHIKLMNPQRSTVWY'*6
            args=SimpleNamespace(threads=1,identification_evalue=1e-5,identification_identity=.3,identification_coverage=.5,competing_ratio=.95)
            for label,headers,expected in [('symbol',['sp|A|one GN=NKX2-3','sp|B|two GN=Nkx2-3'],'homology_assignment'),
                                            ('accession',['sp|A|one','sp|a|two'],'ambiguous_family_assignment')]:
                folder=root/label;folder.mkdir();write_fasta(folder/'refs.faa',[(h,seq) for h in headers])
                prepare(SimpleNamespace(out=str(folder/'db'),fasta=str(folder/'refs.faa')));db,info=verify_db(folder/'db')
                identified,hits=identify(Runner(folder),folder/'identification',{'q':seq},db,info,args)
                self.assertEqual(identified['q']['status'],expected);self.assertEqual(len(identified['q']['near_best']),2)
                if label=='symbol':self.assertEqual(identified['q']['symbol_variants'],['NKX2-3','Nkx2-3'])

class NeighborProposals(unittest.TestCase):
    def test_isoforms_do_not_make_independent_anchors_and_collisions_stay_ambiguous(self):
        from genecontext.stage2 import neighbor_pairs
        models=[dict(model_id=i,genome_id=g,gene_group=k,overlaps_anchor=False) for i,g,k in
                [('a1','a','A'),('a2','a','A'),('b','a','B'),('x','b','X')]]
        def hit(q,s,score):return dict(query=q,subject=s,bitscore=score,identity=.8,coverage_query=1,coverage_subject=1,evalue=1e-20)
        hits=[hit('a1','x',100),hit('a2','x',99),hit('b','x',100),hit('x','a1',100),hit('x','a2',99),hit('x','b',100)]
        args=SimpleNamespace(identification_coverage=.5,identification_identity=.3,competing_ratio=.95)
        pairs=neighbor_pairs(models,hits,args)
        self.assertEqual(len(pairs),4)
        self.assertTrue(all(p['status']=='ambiguous_proposal' for p in pairs))
        self.assertTrue(all(p['kind']=='neighbor_proposal' for p in pairs))
        models[0]['overlaps_anchor']=True
        pairs=neighbor_pairs(models,hits,args)
        self.assertTrue(all(p['kind']=='target_similarity' for p in pairs if 'A' in [p['gene_a'],p['gene_b']]))
    def test_single_annotation_group_reuses_group_without_overwriting_name(self):
        from genecontext.stage2 import assign_groups
        a=dict(model_id='a',source='annotation',genome_id='g',contig='chr',strand='+',gene_group='A',cds_blocks=[[100,200,0]],issues=[])
        p=dict(a,model_id='p',source='tiberius',gene_group='P',proposed_name='UNRESOLVED')
        assign_groups([a,p]);self.assertEqual(p['gene_group'],'A');self.assertEqual(p['proposed_name'],'UNRESOLVED')
