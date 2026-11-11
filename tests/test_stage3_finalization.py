"""Final assignment is a separate, bounded partition of reviewed native models."""
from copy import deepcopy
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
from genecontext.stage3_finalization import SCHEMA, validate_choices, rows_for, prepare
from genecontext.sequence import read_fasta, write_fasta, write_tsv


def fixture():
    models=[]
    for mid,genome in [('a','A'),('b','A'),('c','B')]:
        models.append(dict(model_id=mid,group='g',genome_id=genome,species='Species '+genome,
            original_model='rna_'+mid,original_gene='gene_'+mid,source='annotation',contig='chr1',
            start=10,end=22,strand='+',protein_aa=3,cds_sha256='hash_'+mid,
            quality=dict(coding_status='complete_predicted',source_agreement='exact_CDS_match')))
    identities={m['model_id']:dict(label='NANOG',display_name='NANOG',status='annotation_label',
        basis='source_annotation',ncbi_gene_id='') for m in models}
    identities['b'].update(label='NANOG-like',display_name='NANOG-like',status='curated_locus_label')
    return models,identities


class FinalAssignment(unittest.TestCase):
    def prepared(self,root):
        models,identities=fixture()
        (root/'source.evidence.json').write_text(json.dumps(dict(models=models)))
        folder=root/'queries'/'g';folder.mkdir(parents=True)
        write_tsv(folder/'protein_pairs.tsv',[dict(sequence_a=m['model_id'],sequence_b='reference_q',identity=score,coverage_a=1,coverage_b=1) for m,score in zip(models,[.6,.99,.8])],['sequence_a','sequence_b','identity','coverage_a','coverage_b'])
        for name,seq in [('selected.proteins.faa','MAA'),('selected.cds.fna','ATGGCTGCTTAA')]:write_fasta(root/name,[(m['model_id'],seq) for m in models])
        refs=dict(g=dict(protein=[dict(id='reference_q',symbol='NANOG',accession='NP_test',species='Reference',aliases=[dict(gene_symbol='NANOG')])]))
        payload=dict(groups=[dict(id='g',name='NANOG')],review_id='r',reviewed_pair_count=0)
        with patch('genecontext.stage3_finalization.model_identities',return_value=identities),patch('genecontext.stage3_finalization.export_identities'):
            result=prepare(root,payload,models,refs)
        return result,models

    def test_named_other_copy_never_becomes_automatic_primary_and_exports_exact_native(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);result,models=self.prepared(root)
            self.assertEqual({r['model_id'] for r in result['rows']['primary']},{'a','c'})
            self.assertEqual([r['model_id'] for r in result['rows']['cases']],['b'])
            self.assertEqual(models[1]['copy_identity']['label'],'NANOG-like')
            self.assertTrue(all(r['decision_origin']=='suggestion' for rows in result['rows'].values() for r in rows))
            for prefix,ids in [('primary_targets',{'a','c'}),('special_cases',{'b'})]:
                records=read_fasta(root/(prefix+'.cds.fna'));self.assertEqual(set(records),ids)
                self.assertEqual(set(records.values()),{'ATGGCTGCTTAA'})

    def test_portable_selection_cannot_select_foreign_bucket_reference_or_stale_set(self):
        with tempfile.TemporaryDirectory() as folder:
            result,models=self.prepared(Path(folder));state=result['state'];spec=result['spec'];key=spec['buckets'][0]['key']
            for field,value in [('primary','c'),('primary','reference_q'),('primary',['a','b']),('reason',''),('reviewed','true')]:
                invalid=deepcopy(state);invalid['choices'][key][field]=value
                with self.subTest(field=field,value=value),self.assertRaises(ValueError):validate_choices(invalid,spec)
            for field in ['review_id','set_id']:
                invalid=deepcopy(state);invalid[field]='wrong'
                with self.assertRaises(ValueError):validate_choices(invalid,spec)
            invalid=deepcopy(state);del invalid['choices'][key]
            with self.assertRaises(ValueError):validate_choices(invalid,spec)
            invalid=deepcopy(state);invalid['cases']={'unknown':dict(reason='x',reviewed=True)}
            with self.assertRaises(ValueError):validate_choices(invalid,spec)

    def test_swapping_or_removing_primary_partitions_only_the_selected_set(self):
        with tempfile.TemporaryDirectory() as folder:
            result,models=self.prepared(Path(folder));state=result['state'];spec=result['spec'];key=spec['buckets'][0]['key']
            state['choices'][key]=dict(primary='b',reason='Manual locus-based choice',reviewed=True)
            primary,cases=rows_for(spec,models,validate_choices(state,spec))
            self.assertEqual({r['model_id'] for r in primary},{'b','c'})
            self.assertEqual(next(r['decision_origin'] for r in primary if r['model_id']=='b'),'reviewed')
            state['choices'][key]['primary']='';primary,cases=rows_for(spec,models,state)
            self.assertEqual({r['model_id'] for r in primary},{'c'})
            self.assertEqual({r['model_id'] for r in cases},{'a','b'})
            self.assertEqual(len(primary+cases),3)
            self.assertEqual({r['model_id'] for r in primary+cases},{m['model_id'] for m in models})
