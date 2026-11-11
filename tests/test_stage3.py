"""New STOP-2 scientific decisions and exact predictor-free comparison handoff."""
import copy
import hashlib
import json
from pathlib import Path
import random
import shutil
import tempfile
import unittest

from genecontext.cli import parser
from genecontext.sequence import Model, coding_sequence, write_fasta, sha256
from genecontext.stage2 import model_record, write_outputs
from genecontext.stage2_choices import constraints, validate_choices
from genecontext.stage3 import prepare_inputs, run


def fixture():
    rng=random.Random(39)
    flank=''.join(rng.choice('ACGT') for _ in range(1000))
    dna='ATG'+'GCTGATTGTGAATTTGGTCATATTAAACTGATGAATCCTCAACGTTCTACTGTTTGG'+'TAA'
    models=[];regions=[];sequences={};coding={};proteins={}
    for gid in ['Species_a_GCF_1.1','Species_b_GCF_2.1']:
        rid='region_'+gid
        r=dict(region_id=rid,locus_id=rid,genome_id=gid,gene='TEST',contig='chr',strand='+',
               anchor_start=100,anchor_end=300,window_start=0,window_end=1000,boundary_source='test')
        regions.append(r)
        cds=dna if gid.startswith('Species_a') else dna[:3]+'GCC'+dna[6:]
        sequence=flank[:120]+cds+flank[120+len(cds):];sequences[rid]=sequence
        m=Model('model_'+gid,'chr','+',cds=[(120,120+len(cds),0)],gene='gene_'+gid)
        class Genome:
            def fetch(self,c,a,b,strand='+'):return sequence[a:b]
        nt,aa,status,issues=coding_sequence(m,Genome())
        rec=model_record(m,r,'annotation',nt,aa,status,issues,'TEST')
        rec['identification']=dict(proposed_name='TEST',status='homology_assignment',near_best=[dict(gene_symbol='TEST')])
        rec['target_support']=[dict(subject='q',bitscore=100)]
        models.append(rec);coding[rec['model_id']]=nt;proteins[rec['model_id']]=aa
    data=dict(schema='stage2_review_v1',version='test',review_id='review',source_review_id='parent',
              status='awaiting_manual_review',regions=regions,models=models,neighbor_pairs=[],errors=[],
              query_members=dict(q=[dict(gene_symbol='TEST')]),provenance={})
    state=dict(schema='stage2_draft_v1',review_id='review',regions={r['region_id']:dict(decision='inspect',saved=True,notes='Software control') for r in regions},models={m['model_id']:dict(decision='inspect',role='target',notes='Software control') for m in models},pairs={})
    return data,state,sequences,coding,proteins


def add_paralog(data,state,windows,cds,proteins):
    a=data['models'][0];mid='paralog_copy2';start=420;nt=cds[a['model_id']]
    m=dict(a,model_id=mid,gene_group='copy2_gene',original_gene='copy2_gene',
           original_model='copy2_transcript',start=start,end=start+len(nt),
           cds_blocks=[[start,start+len(nt),0]],overlaps_anchor=False)
    data['models'].append(m);rid=m['region_id'];sequence=windows[rid]
    windows[rid]=sequence[:start]+nt+sequence[start+len(nt):]
    cds[mid]=nt;proteins[mid]=proteins[a['model_id']]
    state['models'][mid]=dict(decision='inspect',role='target',selection_kind='additional_copy',
        comparison_group='TEST',copy_label='TEST_copy2',notes='Explicitly selected additional copy')
    return m


class NewReviewChoices(unittest.TestCase):
    def test_added_paralog_has_own_anchor_but_query_and_family_checks_still_apply(self):
        d,s,windows,cds,aa=fixture();m=add_paralog(d,s,windows,cds,aa)
        self.assertEqual(validate_choices(d,s),[])
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name,records in [('windows.fna',windows),('cds.fna',cds),('proteins.faa',aa)]:write_fasta(root/name,records.items())
            selected=prepare_inputs(d,s,root);self.assertEqual(len(selected),3)
            row=next(r for r,_,_,_ in selected if r['model_id']==m['model_id'])
            self.assertEqual((row['target_start'],row['target_end']),(m['start'],m['end']))
            self.assertEqual(row['copy_label'],'TEST_copy2')
        m['target_support']=[]
        self.assertTrue(any('query support' in x['message'] for x in validate_choices(d,s)))
        m['target_support']=[dict(subject='q',bitscore=100)]
        m['identification']=dict(near_best=[dict(gene_symbol='OTHER')])
        self.assertTrue(any('another gene' in x['message'] for x in validate_choices(d,s)))

    def test_valid_models_and_explicit_override(self):
        d,s,*_=fixture();self.assertEqual(validate_choices(d,s),[])
        m=d['models'][0];m['overlaps_anchor']=False;m['target_support']=[]
        self.assertTrue(any('override' in x['message'] for x in validate_choices(d,s)))
        c=s['models'][m['model_id']];c['target_override']=True
        self.assertEqual(validate_choices(d,s),[])
        c['target_override']='true';self.assertTrue(validate_choices(d,s))

    def test_conflicting_identity_and_malformed_records_are_rejected(self):
        d,s,*_=fixture();d['models'][0]['identification']['near_best']=[dict(gene_symbol='OTHER')]
        self.assertTrue(any('another gene' in x['message'] for x in validate_choices(d,s)))
        for value in [None, [], {'schema':'wrong'}]:self.assertTrue(validate_choices(d,value))
        s['models'][d['models'][0]['model_id']]=None
        self.assertTrue(validate_choices(d,s))

    def test_duplicate_physical_gene_and_comparison_copies(self):
        d,s,*_=fixture();a=d['models'][0];other=dict(a,model_id='repeat')
        d['models'].append(other);s['models']['repeat']=dict(s['models'][a['model_id']])
        self.assertTrue(any('more than once' in x['message'] for x in validate_choices(d,s)))
        other.update(gene_group='distinct',original_gene='other_gene',start=400,end=500,cds_blocks=[[400,500,0]])
        self.assertEqual(validate_choices(d,s),[])
        s['models'][a['model_id']]['copy_label']='copy1';s['models']['repeat']['copy_label']='copy1'
        self.assertTrue(any('distinct copy labels' in x['message'] for x in validate_choices(d,s)))
        s['models']['repeat']['comparison_group']='TEST_copy2'
        self.assertEqual(validate_choices(d,s),[])

    def test_distinct_annotation_ids_do_not_merge_and_old_window_variants_do(self):
        d,s,*_=fixture();a=d['models'][0]
        duplicate=dict(a,model_id='alternative',gene_group='other_group',region_id='another_window',source='tiberius',cds_blocks=[[120,135,0],[140,150,0],[160,a['end'],0]])
        # Identical outer bounds + substantial CDS overlap reconcile an old
        # annotation-free observation, without relying on the window ID.
        a['source']='tiberius';d['models'].append(duplicate)
        meta=constraints(d);self.assertEqual(meta[a['model_id']]['physical_id'],meta['alternative']['physical_id'])
        a['source']='annotation';duplicate['source']='annotation';duplicate['original_gene']='different'
        meta=constraints(d);self.assertNotEqual(meta[a['model_id']]['physical_id'],meta['alternative']['physical_id'])

    def test_neighbor_pair_cannot_survive_excluded_regions(self):
        d,s,*_=fixture();a,b=d['models'];p=dict(gene_a=a['gene_group'],gene_b=b['gene_group'],model_a=a['model_id'],model_b=b['model_id'],kind='neighbor_proposal',decision='inspect')
        d['neighbor_pairs']=[p];s['pairs'][':'.join(sorted([p['gene_a'],p['gene_b']]))]=p
        for c in s['regions'].values():c['decision']='exclude'
        s['models']={}
        self.assertTrue(any('neighbor pair' in x['message'] for x in validate_choices(d,s)))

    def test_neighbor_model_alternatives_require_one_consistent_choice(self):
        d,s,*_=fixture();a=d['models'][0]
        for name,digest in [('n1','one'),('n2','two')]:
            m=dict(a,model_id=name,gene_group='neighbor',original_gene='neighbor',
                   start=600,end=700,cds_blocks=[[600,700,0]],cds_sha256=digest,overlaps_anchor=False)
            d['models'].append(m);s['models'][name]=dict(decision='inspect',role='neighbor',notes='Software control')
        self.assertTrue(any('Conflicting coding models' in x['message'] for x in validate_choices(d,s)))
        del s['models']['n2'];self.assertEqual(validate_choices(d,s),[])

    def test_prepare_inputs_exact_native_roundtrip_and_tamper_detection(self):
        d,s,windows,cds,aa=fixture()
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name,records in [('windows.fna',windows),('cds.fna',cds),('proteins.faa',aa)]:write_fasta(root/name,records.items())
            prepared=prepare_inputs(d,s,root);self.assertEqual(len(prepared),2)
            mid=d['models'][0]['model_id'];aa[mid]='M'+aa[mid]
            write_fasta(root/'proteins.faa',aa.items())
            with self.assertRaisesRegex(ValueError,'hash or integrity'):prepare_inputs(d,s,root)


@unittest.skipUnless(shutil.which('mafft') and shutil.which('blastn'),'requires real MAFFT and BLASTN')
class RealStage3(unittest.TestCase):
    def test_reviewed_models_reach_actual_alignments_and_catalogue_without_predictor(self):
        d,s,windows,cds,aa=fixture();add_paralog(d,s,windows,cds,aa)
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'source';source.mkdir()
            for name,records in [('windows.fna',windows),('cds.fna',cds),('proteins.faa',aa)]:write_fasta(source/name,records.items())
            # Stage-2 BLAST queries are all extracted proteins, while subjects
            # are the supplied target panel. Only the latter may guide MSA.
            (source/'target_support').mkdir()
            write_fasta(source/'target_support/queries.faa',{'unrelated_neighbor':'MACDEFGHIKLMNPQRSTVWY'}.items())
            write_fasta(source/'target_support/subjects.faa',{'q':next(iter(aa.values()))}.items())
            (source/'evidence.json').write_text(json.dumps(d));write_outputs(source,d)
            decisions=root/'choices.json';decisions.write_text(json.dumps(s))
            args=parser().parse_args(['stage3','--results',str(source),'--decisions',str(decisions),'--out',str(root/'result'),'--threads','1','--control-run'])
            args.validate_only=True;self.assertEqual(run(args),0);self.assertFalse((root/'result').exists())
            args.validate_only=False;self.assertEqual(run(args),0)
            manifest=json.loads((root/'result/manifest.json').read_text());self.assertTrue(manifest['control_run'])
            self.assertEqual(manifest['selected_models'],3)
            self.assertEqual(manifest['selected_species'],2)
            commands=[Path(c['argv'][0]).name for c in manifest['commands']]
            self.assertIn('mafft',commands);self.assertIn('blastn',commands)
            self.assertFalse(any(x in commands for x in ['tblastn','blastp','miniprot','tiberius']))
            self.assertTrue(list((root/'result/queries').glob('*/cds.codon.aligned.fna')))
            from genecontext.sequence import read_fasta
            aligned=read_fasta(next((root/'result/queries').glob('*/proteins.aligned.faa')))
            self.assertEqual(set(aligned),set(aa)|{'reference_q'})
            profile=next((root/'result/queries').glob('*/protein_profile.tsv'))
            import csv
            with profile.open() as handle:
                self.assertEqual({r['eligible_loci'] for r in csv.DictReader(handle,delimiter='\t')},{'3'})
            self.assertIn('TEST_copy2',(root/'result/report.html').read_text())
            self.assertIn('SOFTWARE CONTROL',(root/'result/report.html').read_text())
            self.assertTrue((root/'result/context.html').exists())
            self.assertTrue((root/'result/selected.upstream.fna').exists())
            self.assertTrue((root/'result/model_quality.tsv').exists())
            with (root/'result/species_conservation_profile.tsv').open() as handle:
                self.assertEqual({r['eligible_species'] for r in csv.DictReader(handle,delimiter='\t')},{'2'})
            from genecontext.stage3_report import refresh
            args_report=parser().parse_args(['stage3-report','--results',str(root/'result'),'--out',str(root/'refreshed')])
            self.assertEqual(refresh(args_report),0)
            for file in ['selected.cds.fna','selected.proteins.faa','decisions.snapshot.json','comparison.json','manifest.json']:
                self.assertEqual(sha256(root/'result'/file),sha256(root/'refreshed'/file))
            # Source evidence and review identity remain sealed and unchanged.
            self.assertTrue(all(sha256(source/n)==v for n,v in json.loads((source/'seal.json').read_text()).items()))
            wrong=copy.deepcopy(s);wrong['review_id']='another';decisions.write_text(json.dumps(wrong));args.out=str(root/'wrong')
            with self.assertRaisesRegex(ValueError,'another stage 2'):run(args)
            self.assertFalse((root/'wrong').exists())
