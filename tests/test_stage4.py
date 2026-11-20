"""STOP 4 must preserve reviewed sequence identity and expose missing evidence."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.cli import parser
from genecontext.sequence import write_fasta,write_tsv,sha256,translate
from genecontext.stage4 import run,checked_source
from genecontext.stage4_phylogeny import unique_alignment,parse_newick,tip_names,tree_svg,build_tree
from genecontext.stage4_structure import ca_atoms,import_structures,parse_usalign,compare_structures,RESIDUES


def digest(text):return hashlib.sha256(text.encode()).hexdigest()


def fixture(root,reviewed=False):
    proteins={'mA':'MA'+'D'*35,'mB':'MA'+'E'*35}
    cds={k:'ATG'+'GCT'+('GAT' if k=='mA' else 'GAA')*35+'TAA' for k in proteins}
    rows=[];selected=[];q=[];evidence={};buckets=[];choices={}
    for mid in proteins:
        row=dict(target='TEST',species='Species '+mid,genome_id=mid+'_GCF_1.1',model_id=mid,
            assembly_accession='GCF_1.1',physical_id=mid,original_model='rna_'+mid,original_gene='gene_TEST',
            source='annotation',contig='chr1',start=10,end=10+len(cds[mid]),strand='+',
            protein_aa=len(proteins[mid]),protein_sha256=digest(proteins[mid]),cds_sha256=digest(cds[mid]),
            coding_status='complete_predicted',window_edge=False,protein_file='sequences/'+mid+'.protein.faa',cds_file='sequences/'+mid+'.cds.fna')
        rows.append(row);selected.append(dict(model_id=mid,group='g',genome_id=row['genome_id']))
        quality=dict(model_id=mid,coding_status='complete_predicted',source_agreement='different_CDS_models',cautions='annotation_prediction_disagreement')
        q.append(quality)
        evidence[mid]=dict(copy_identity=dict(display_name='TEST',label='TEST',basis='provided_annotation',status='annotation_label',ncbi_gene_id=''),assignment=dict(best_target=dict(reference_id='reference_q',accession='NP_1',symbol='TEST',identity=.8,coverage_model=1.,coverage_reference=1.),known_other_copy=False,alternate_preferred=False,complete=True,facts=['Fixture']))
        buckets.append(dict(key=mid,gene='TEST',group='g',genome_id=row['genome_id'],species=row['species'],model_ids=[mid]))
        choices[mid]=dict(primary=mid,reason='Explicit test choice',reviewed=reviewed)
    spec=dict(review_id='r',set_id=digest(json.dumps(['r',sorted([('g',r['genome_id'],r['model_id'],r['cds_sha256']) for r in rows])],sort_keys=True)),buckets=buckets)
    state=dict(schema='stage3_final_choices_v1',review_id='r',set_id=spec['set_id'],choices=choices,cases={})
    for name,value in [('manifest.json',dict(stage='stage3_comparison',status='completed',source_review_id='r',errors=[])),('comparison.json',dict(groups=[dict(id='g',name='TEST',queries=['q'])],selected=selected)),('finalization.evidence.json',dict(spec=spec,models=evidence)),('finalization.decisions.json',state)]:
        (root/name).write_text(json.dumps(value))
    for name,value in [('sequence_catalogue.tsv',rows),('model_quality.tsv',q)]:write_tsv(root/name,value,list(value[0]))
    refs=[dict(model_id=mid,kind='protein',reference_id='reference_q',reference_accession='NP_1',reference_species='Reference species',identity=.8,coverage_model=1.,coverage_reference=1.) for mid in proteins]
    write_tsv(root/'reference_comparisons.tsv',refs,list(refs[0]))
    write_tsv(root/'synteny.tsv',[],['locus_a','locus_b','unique_anchors'])
    write_fasta(root/'selected.proteins.faa',proteins.items());write_fasta(root/'selected.cds.fna',cds.items())
    folder=root/'queries/g';folder.mkdir(parents=True)
    write_fasta(folder/'proteins.aligned.faa',list(proteins.items())+[('reference_q',proteins['mA'])])
    seal(root)
    return state,proteins


def seal(root):
    (root/'seal.json').write_text(json.dumps({str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'}))


def pdb(root,sequence):
    reverse={v:k for k,v in RESIDUES.items()};path=root/'synthetic.pdb'
    path.write_text('\n'.join(f'ATOM  {n:5d}  CA  {reverse[aa]:3s} A{n:4d}    {n*1.0:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{85.:6.2f}           C' for n,aa in enumerate(sequence,1))+'\nEND\n')
    return path


class Stage4Contract(unittest.TestCase):
    def args(self,source,out,*extra):
        return parser().parse_args(['stage4','--results',str(source),'--out',str(out),'--tree','none',*extra])

    def test_preview_is_explicit_no_orthology_or_mrna_claim_and_exact_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();_,proteins=fixture(source);out=Path(tmp)/'out'
            with self.assertRaisesRegex(ValueError,'not reviewed'):run(self.args(source,out))
            self.assertFalse(out.exists())
            run(self.args(source,out,'--preview'))
            data=json.loads((out/'evidence.json').read_text())
            self.assertTrue(data['preview']);self.assertEqual(data['counts']['primary'],2)
            self.assertTrue(all(r['model_assessment']=='Check model' and r['synteny_status']=='Not resolved' and r['structure'] is None for r in data['rows']))
            self.assertTrue(all('orthology' in ' '.join(r['limitations']) for r in data['rows']))
            from genecontext.sequence import read_fasta
            self.assertEqual(read_fasta(out/'primary_targets.proteins.faa'),proteins)
            self.assertFalse(json.loads((out/'manifest.json').read_text())['biological_approval_asserted'])
            checked_source(source) # Never edit input seal or decisions.

    def test_bad_seal_foreign_choice_and_protein_cds_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();state,_=fixture(source,True);out=Path(tmp)/'out'
            (source/'selected.proteins.faa').write_text('>mA\nMA\n')
            with self.assertRaisesRegex(ValueError,'seal mismatch'):run(self.args(source,out))
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();state,_=fixture(source,True)
            state['choices']['mA']['primary']='mB';(source/'finalization.decisions.json').write_text(json.dumps(state));seal(source)
            with self.assertRaisesRegex(ValueError,'own gene/genome bucket'):run(self.args(source,Path(tmp)/'out'))
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();fixture(source,True)
            p=source/'selected.cds.fna';p.write_text(p.read_text().replace('GAT','GAA'))
            # Even a consistently rehashed CDS must still translate to the chosen protein.
            from genecontext.sequence import read_fasta
            coding=read_fasta(p)
            rows=__import__('genecontext.stage3_analysis',fromlist=['read_table']).read_table(source/'sequence_catalogue.tsv')
            for r in rows:r['cds_sha256']=digest(coding[r['model_id']])
            write_tsv(source/'sequence_catalogue.tsv',rows,list(rows[0]));seal(source)
            with self.assertRaisesRegex(ValueError,'translation differs'):run(self.args(source,Path(tmp)/'out'))

    def test_stored_msa_cannot_substitute_a_different_selected_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();fixture(source,True)
            p=source/'queries/g/proteins.aligned.faa';p.write_text(p.read_text().replace('D','E'));seal(source)
            with self.assertRaisesRegex(ValueError,'Stored MSA differs'):run(self.args(source,Path(tmp)/'out'))

    def test_reviewed_primary_and_separate_cases_remain_partitioned(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'source';source.mkdir();state,_=fixture(source,True)
            state['choices']['mA']['primary']='';(source/'finalization.decisions.json').write_text(json.dumps(state));seal(source)
            out=Path(tmp)/'out';run(self.args(source,out))
            data=json.loads((out/'evidence.json').read_text());self.assertFalse(data['preview'])
            self.assertEqual(data['counts']['primary'],1);self.assertEqual(data['counts']['special_cases'],1)
            self.assertEqual({r['model_id'] for r in data['rows']},{'mA','mB'})

    def test_tree_dedup_mask_aliases_and_foreign_root(self):
        aln={'a':'MA---DD','b':'MA---DD','c':'MEAAAED','d':'MDAAADE'}
        masked,leaves,mask=unique_alignment(aln,{k:k for k in aln},.75)
        self.assertEqual(len(leaves),3);self.assertEqual(leaves[0]['sequence_ids'],['a','b'])
        self.assertEqual(mask['retained_columns_1based'],[1,2,6,7])
        self.assertEqual(len(masked),3)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'both ingroup and outgroup'):
                build_tree(Path(tmp),aln,{k:k for k in aln},None,outgroups=['a'],root_basis='Test')
        svg=tree_svg('(T00001:0.1,(T00002:0.2,T00003:0.3)98:0.1);',leaves)
        self.assertIn('identical',svg)
        with self.assertRaises(ValueError):tree_svg('(T00001:0.1,foreign:0.2);',leaves)
        self.assertEqual(tip_names(parse_newick('(a:.1,b:0);')),['a','b'])

    def test_structures_require_exact_sequence_and_declared_plddt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq='M'+'A'*25;path=pdb(root,seq)
            self.assertEqual(''.join(a['aa'] for a in ca_atoms(path)),seq)
            entry=dict(sequence_id='m',source='Synthetic test only',file=path.name,file_sha256=sha256(path),protein_sha256=digest(seq))
            cat=root/'catalogue.json';cat.write_text(json.dumps(dict(schema='stage4_structures_v1',entries=[entry])))
            imported,_=import_structures(cat,{'m':seq},root);self.assertIsNone(imported['m']['mean_plddt'])
            entry['confidence']='plddt_bfactor';cat.write_text(json.dumps(dict(schema='stage4_structures_v1',entries=[entry])))
            imported,_=import_structures(cat,{'m':seq},root);self.assertEqual(imported['m']['fraction_plddt_ge70'],1.)
            entry['protein_sha256']=digest(seq+'A');cat.write_text(json.dumps(dict(schema='stage4_structures_v1',entries=[entry])))
            with self.assertRaisesRegex(ValueError,'CA sequence'):import_structures(cat,{'m':seq+'A'},root)

    def test_usalign_both_normalizations_and_no_fabricated_missing_metrics(self):
        text='Aligned length= 29, RMSD= 1.25, Seq_ID=n_identical/n_aligned= 0.9\nTM-score= 0.87 (normalized by length of Structure_1: L=30)\nTM-score= 0.81 (normalized by length of Structure_2: L=35)\n'
        value=parse_usalign(text);self.assertEqual((value['tm_score_a'],value['tm_score_b']),(.87,.81))
        with self.assertRaises(ValueError):parse_usalign(text.splitlines()[0])
        with tempfile.TemporaryDirectory() as tmp:
            rows=compare_structures({}, {},Path(tmp),None)
            self.assertEqual(rows,[])


if __name__=='__main__':unittest.main()
