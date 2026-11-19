"""Offline protocol controls; ordinary tests never contact AFDB or spend credits."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs,urlsplit

from genecontext.alphafold_db import BASE,Client,collect
from genecontext.af3_requests import prepare
from genecontext.cli import parser
from genecontext.sequence import sha256
from genecontext.stage4 import run
from genecontext.stage4_structure import import_structures
from test_stage4 import fixture,pdb


def hashes(root):
    return {str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file()}


class Protocol:
    def __init__(self,root,sequences):
        self.calls=[];self.routes={};self.sequences={};self.metadata={}
        for i,seq in enumerate(sequences,1):
            mid='AF-TEST'+str(i)+'-F1';digest=hashlib.md5(seq.encode()).hexdigest()
            self.sequences[digest]=seq
            summary=dict(entry=dict(sequence=seq),structures=[dict(summary=dict(model_identifier=mid,
                provider='AlphaFold DB',coverage=1,sequence_identity=1,oligomeric_state='MONOMER',created='2026-01-01',confidence_avg_local_score=85.))])
            self.routes[digest]=summary
            self.metadata[mid]=[dict(modelEntityId=mid,providerId='GDM',isComplex=False,chainId='A',sequence=seq,
                sequenceStart=1,sequenceEnd=len(seq),latestVersion=6,uniprotAccession='TEST'+str(i),
                pdbUrl=BASE+'/files/'+mid+'.pdb',paeDocUrl=BASE+'/files/'+mid+'.pae.json')]
            self.routes[mid+'.pdb']=pdb(root,seq).read_bytes()
            self.routes[mid+'.pae.json']=[dict(predicted_aligned_error=[[1.]*len(seq) for _ in seq])]

    def __call__(self,url,maximum):
        self.calls.append(url);parsed=urlsplit(url)
        if parsed.path.endswith('/sequence/summary'):
            key=parse_qs(parsed.query)['id'][0];data=self.routes.get(key)
        elif '/prediction/' in parsed.path:
            data=self.metadata[parsed.path.rsplit('/',1)[1]]
        else:data=self.routes[parsed.path.rsplit('/',1)[1]]
        return data if isinstance(data,bytes) or data is None else json.dumps(data).encode()


class ExactAFDB(unittest.TestCase):
    def test_exact_sequence_coordinates_pae_dedup_and_offline_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq='M'+'A'*25;protocol=Protocol(root,[seq])
            catalogue,rows=collect({'model_a':seq,'same_protein_query':seq},root/'fetch',root/'cache',transport=protocol)
            self.assertEqual(len(protocol.calls),4)
            self.assertTrue(all(r['status']=='available' for r in rows))
            data=json.loads(catalogue.read_text());self.assertEqual(len(data['entries']),2)
            self.assertEqual(data['entries'][0]['file'],data['entries'][1]['file'])
            destination=root/'imported';destination.mkdir()
            imported,_=import_structures(catalogue,{'model_a':seq,'same_protein_query':seq},destination)
            self.assertEqual(imported['model_a']['mean_plddt'],85)
            self.assertEqual(imported['model_a']['pae_status'],'available_not_interpreted')
            _,rows=collect({'model_a':seq},root/'offline',root/'cache',offline=True,transport=lambda *_:self.fail('Offline request'))
            self.assertEqual(rows[0]['status'],'available')
            cache_pdb=next((root/'cache').glob('*/model.pdb'));cache_pdb.write_text('corrupt')
            with self.assertRaisesRegex(ValueError,'changed|hash|checksum|Cache'):
                Client(root/'cache',offline=True).fetch(seq)

    def test_sequence_or_coordinate_mismatch_never_imports_a_proxy(self):
        for mismatch in ['metadata','coordinates']:
            with self.subTest(mismatch=mismatch),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);seq='M'+'A'*25;protocol=Protocol(root,[seq])
                if mismatch=='metadata':protocol.metadata['AF-TEST1-F1'][0]['sequence']='M'+'G'*25
                else:protocol.routes['AF-TEST1-F1.pdb']=pdb(root,'M'+'G'*25).read_bytes()
                catalogue,rows=collect({'m':seq},root/'fetch',root/'cache',transport=protocol)
                self.assertEqual(json.loads(catalogue.read_text())['entries'],[])
                self.assertEqual(rows[0]['status'],'not_found')

    def test_optional_bad_pae_keeps_verified_coordinates_without_inventing_pae(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq='M'+'A'*25;protocol=Protocol(root,[seq]);protocol.routes['AF-TEST1-F1.pae.json']=[dict(predicted_aligned_error=[[1.]])]
            catalogue,rows=collect({'m':seq},root/'fetch',root/'cache',transport=protocol)
            self.assertEqual(rows[0]['status'],'available');self.assertEqual(rows[0]['pae_status'],'unavailable')
            self.assertNotIn('pae_file',json.loads(catalogue.read_text())['entries'][0])

    def test_untrusted_download_urls_rejected_and_unsupported_or_limited_inputs_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq='M'+'A'*25;protocol=Protocol(root,[seq]);protocol.metadata['AF-TEST1-F1'][0]['pdbUrl']='http://127.0.0.1/private'
            _,rows=collect({'m':seq,'second':'MG','bad':'MX'},root/'fetch',root/'cache',limit=1,transport=protocol)
            self.assertEqual(rows[0]['status'],'unavailable');self.assertIn('URL',rows[0]['reason'])
            self.assertEqual(rows[1]['status'],'not_requested')
            self.assertFalse(any('127.0.0.1' in url for url in protocol.calls))
            with self.assertRaisesRegex(ValueError,'canonical'):Client(root/'cache').fetch('MX')

    def test_offline_cache_miss_never_submits_or_downloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            catalogue,rows=collect({'m':'MA'},root/'fetch',root/'cache',offline=True,transport=lambda *_:self.fail('Network'))
            self.assertEqual(rows[0]['status'],'offline_cache_miss')
            self.assertFalse(json.loads((catalogue.parent/'retrieval.json').read_text())['prediction_submitted'])

    def test_malformed_service_response_is_explicit_unavailability(self):
        for response in [[],dict(entry=None,structures=[]),dict(entry=dict(sequence='MA'),structures=[None])]:
            with self.subTest(response=response),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                _,rows=collect({'m':'MA'},root/'fetch',root/'cache',transport=lambda *_:json.dumps(response).encode())
                self.assertEqual(rows[0]['status'],'unavailable')

    def test_native_stop4_opt_in_uses_importer_and_keeps_source_and_choices_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();_,proteins=fixture(source,True);before=hashes(source)
            protocol=Protocol(root,list(proteins.values()));out=root/'stop4'
            args=parser().parse_args(['stage4','--results',str(source),'--out',str(out),'--tree','none',
                '--structure-provider','afdb','--structure-cache',str(root/'cache')])
            with patch('genecontext.alphafold_db.Download',return_value=protocol):run(args)
            evidence=json.loads((out/'evidence.json').read_text())
            self.assertTrue(all(r['structure']['exact_sequence_match'] for r in evidence['rows']))
            self.assertEqual(evidence['structure_viewer']['unique_structures'],2)
            self.assertEqual(hashes(source),before)
            self.assertTrue((out/'structure_viewer.html').exists())

    def test_native_stop4_default_never_contacts_a_prediction_or_database_service(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();fixture(source,True)
            args=parser().parse_args(['stage4','--results',str(source),'--out',str(root/'stop4'),'--tree','none'])
            with patch('genecontext.alphafold_db.urlopen',side_effect=AssertionError('No network')):run(args)


class AF3Requests(unittest.TestCase):
    def test_cloud_compatible_jobs_bind_full_sequence_and_never_submit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with patch('genecontext.alphafold_db.urlopen',side_effect=AssertionError('No network')):
                rows=prepare({'selected':'MA','same_sequence':'MA','later':'MG'},root/'prepared',limit=1)
            self.assertEqual(rows[0]['status'],'awaiting_prediction');self.assertEqual(rows[1]['input_file'],rows[0]['input_file'])
            self.assertEqual(rows[2]['status'],'not_requested')
            data=json.loads((root/'prepared'/rows[0]['input_file']).read_text())
            self.assertEqual(data['sequences'],[dict(protein=dict(id='A',sequence='MA'))])
            self.assertEqual((data['dialect'],data['version'],data['modelSeeds']),('alphafold3',4,[42]))
            server=json.loads((root/'prepared'/rows[0]['server_input_file']).read_text())
            self.assertEqual(server[0]['sequences'],[dict(proteinChain=dict(sequence='MA',count=1))])
            request=json.loads((root/'prepared/requests.json').read_text());self.assertFalse(request['prediction_submitted'])
            self.assertEqual(request['requests'][0]['source_sequence_ids'],['selected','same_sequence'])


if __name__=='__main__':unittest.main()
