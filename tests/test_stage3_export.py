"""Bulk default exports preserve native models, literal comments and review state."""
from copy import deepcopy
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from genecontext.sequence import sha256
from genecontext.stage3_export import write_default_export,write_collection,query_export


def fixture(root):
    source=root/'stage2';source.mkdir()
    data=dict(source_review_id='stop1',regions=[dict(region_id='r',genome_id='genome',
        reason='other',notes='Стоп 1\nКомментарий\tс табом')],models=[dict(model_id=mid,region_id='r') for mid in ['a','b']])
    content=json.dumps(data,ensure_ascii=False)
    (source/'evidence.json').write_text(content)
    draft=dict(review_id='stop1',pairs=[dict(genome_id='genome',notes='Исходный комментарий\nцеликом',selected=['r'])])
    (source/'selection.input.json').write_bytes((json.dumps(draft,ensure_ascii=False,indent=2)+'\n').replace('\n','\r\n').encode())
    (source/'seal.json').write_text(json.dumps({'selection.input.json':sha256(source/'selection.input.json')}))
    report=root/'report';report.mkdir();(report/'source.evidence.json').write_text(content)
    state=dict(review_id='review',models={'a':dict(notes='Модель A'),'b':dict(notes='Копия B')},
        regions={'r':dict(notes='Решение STOP2')},automatic_comments={'r':['Automatic warning']},
        pairs={'p':dict(notes='Pair comment retained')},review_plan={'r':dict(reason='Plan note')})
    (report/'decisions.snapshot.json').write_text(json.dumps(state,ensure_ascii=False))
    payload=dict(review_id='review',models=[dict(model_id=mid) for mid in ['a','b']],parameters=dict(results=str(source)))
    choice=dict(schema='stage3_final_choices_v1',review_id='review',set_id='exact_set',choices={'bucket':dict(primary='a',reason='Default proposal',reviewed=False)},cases={})
    rows=[dict(gene='GENE',model_id=mid,classification=role,selection_reason=reason,decision_origin='suggestion') for mid,role,reason in [('a','primary_target','Default proposal'),('b','special_case','Additional copy')]]
    final=dict(state=choice,rows=dict(primary=rows[:1],cases=rows[1:]),raw_sequences=dict(protein=dict(a='MAA',b='MCC'),cds=dict(a='ATGGCTGCTTAA',b='ATGTGTTGTTAA')))
    return source,report,payload,final


class DefaultExports(unittest.TestCase):
    def test_queries_are_exported_exactly_with_source_identity_and_do_not_enter_choices(self):
        with tempfile.TemporaryDirectory() as directory:
            source,report,payload,final=fixture(Path(directory))
            folder=report/'queries/g';folder.mkdir(parents=True)
            # Identical sequences with different IDs must both survive.
            (folder/'proteins.input.faa').write_text('>a\nMAA\n>reference_q1\nMCC\n>reference_q2\nMCC\n')
            payload.update(groups=[dict(id='g',name='GENE',queries=['q1','q2'])],
                query_members={'q1':[dict(protein='NP_test.1',species='Species Q',gene_symbol='GENE')]},
                alignments={'g':{'protein':{'reference_q1':'M-CC','reference_q2':'MCC-'}}})
            before=deepcopy(final);result=write_default_export(report,payload,final)
            self.assertEqual(final,before)
            self.assertEqual(result['record']['query_proteins'],2)
            with zipfile.ZipFile(report/result['record']['package']) as archive:
                query=archive.read('query_references.faa').decode()
                self.assertEqual(query,'>reference_q1\nMCC\n>reference_q2\nMCC\n')
                combined=archive.read('all_selected_with_queries.proteins.faa').decode()
                self.assertEqual(combined.count('>'),4)
                self.assertEqual(archive.read('primary_targets_with_queries.proteins.faa').decode().count('>'),3)
                refs=list(csv.DictReader(io.StringIO(archive.read('query_references.tsv').decode()),delimiter='\t'))
                self.assertEqual(refs[0]['accession'],'NP_test.1')
                self.assertEqual(refs[0]['source_file_sha256'],sha256(folder/'proteins.input.faa'))
                self.assertEqual(refs[0]['protein_sha256'],hashlib.sha256(b'MCC').hexdigest())
                self.assertEqual(refs[0]['cds_status'],'not_supplied_to_stop3')
                self.assertEqual(json.loads(archive.read('stage3.final.decisions.json')),before['state'])
                self.assertNotIn('query_references.cds.fna',archive.namelist())

    def test_missing_or_conflicting_queries_cannot_silently_disappear(self):
        with tempfile.TemporaryDirectory() as directory:
            source,report,payload,final=fixture(Path(directory))
            payload['groups']=[dict(id='g',name='GENE',queries=['q'])]
            with self.assertRaisesRegex(ValueError,'input is missing'):query_export(report,payload)
            folder=report/'queries/g';folder.mkdir(parents=True)
            (folder/'proteins.input.faa').write_text('>a\nMAA\n')
            with self.assertRaisesRegex(ValueError,'protein is missing'):query_export(report,payload)
            (folder/'proteins.input.faa').write_text('>reference_q\nMAA\n')
            payload['alignments']={'g':{'protein':{'reference_q':'MCC'}}}
            with self.assertRaisesRegex(ValueError,'alignment differs'):query_export(report,payload)
            payload.pop('alignments')
            payload['groups'].append(dict(id='g2',name='GENE',queries=['q']))
            folder=report/'queries/g2';folder.mkdir()
            (folder/'proteins.input.faa').write_text('>reference_q\nMCC\n')
            with self.assertRaisesRegex(ValueError,'Conflicting query'):query_export(report,payload)

    def test_one_package_preserves_defaults_all_models_and_literal_previous_comments(self):
        with tempfile.TemporaryDirectory() as directory:
            source,report,payload,final=fixture(Path(directory));before=deepcopy(final)
            result=write_default_export(report,payload,final)
            self.assertEqual(final,before)
            with zipfile.ZipFile(report/result['record']['package']) as archive:
                state=json.loads(archive.read('stage3.final.decisions.json'))
                self.assertFalse(state['choices']['bucket']['reviewed'])
                self.assertEqual(archive.read('comments/stop1.selection.input.json'),(source/'selection.input.json').read_bytes())
                self.assertIn('Pair comment retained',archive.read('comments/stop2.decisions.json').decode())
                rows=list(csv.DictReader(io.StringIO(archive.read('models_with_comments.tsv').decode()),delimiter='\t'))
                self.assertEqual({r['model_id'] for r in rows},{'a','b'})
                self.assertEqual(rows[0]['stop1_notes'],'Стоп 1\nКомментарий\tс табом')
                self.assertIn('Automatic warning',rows[0]['stop2_automatic_comments'])
                self.assertIn('Копия B',rows[1]['stop2_model_notes'])
                self.assertIn(b'ATGGCTGCTTAA',archive.read('primary_targets.cds.fna'))
                self.assertIn(b'ATGTGTTGTTAA',archive.read('special_cases.cds.fna'))
                self.assertEqual(archive.testzip(),None)

    def test_changed_prior_snapshot_is_rejected_and_missing_history_is_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            source,report,payload,final=fixture(Path(directory))
            (source/'selection.input.json').write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):write_default_export(report,payload,final)
            payload['parameters']={}
            result=write_default_export(report,payload,final)
            self.assertEqual(result['history']['stop1_history_status'],'inherited_region_comments_only')
            # The portable history remains usable when the original server is unavailable.
            payload['parameters']['results']='missing'
            self.assertEqual(write_default_export(report,payload,final)['history'],result['history'])

    def test_collection_namespaces_repeated_model_ids_and_preserves_every_gene_choice(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);reports=[]
            for key in ['GENE1','GENE2']:
                folder=root/key;folder.mkdir();source,report,payload,final=fixture(folder)
                write_default_export(report,payload,final)
                reports.append(dict(key=key,results=str(report),url=key+'/report.html'))
            result=write_collection(root/'collection',reports)
            self.assertEqual((result['reports'],result['models']),(2,4))
            with zipfile.ZipFile(root/'collection/all_default_models_with_comments.zip') as archive:
                for key in ['GENE1','GENE2']:self.assertIn(key+'/stage3.final.decisions.json',archive.namelist())
                bundle=json.loads(archive.read('stage3.defaults.bundle.json'))
                self.assertTrue(all(not x['choices']['choices']['bucket']['reviewed'] for x in bundle['reports']))
                self.assertEqual(len(list(csv.DictReader(io.StringIO(archive.read('all_models_with_comments.tsv').decode()),delimiter='\t'))),4)
            with self.assertRaisesRegex(ValueError,'Duplicate'):write_collection(root/'duplicate',[reports[0],reports[0]])
            (Path(reports[0]['results'])/'default_models_with_comments.zip').write_text('broken')
            with self.assertRaisesRegex(ValueError,'changed'):write_collection(root/'tampered',reports)
