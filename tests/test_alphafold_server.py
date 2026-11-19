"""Offline Server ZIP controls; no browser sessions or predictions in tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from genecontext.af3_requests import prepare
from genecontext.alphafold_server import import_result
from genecontext.sequence import sha256
from genecontext.stage4_structure import import_structures
from test_stage4 import pdb

HAS_GEMMI=importlib.util.find_spec('gemmi') is not None


def fixture(root,sequence='MAAAA'):
    import gemmi
    rows=prepare({'selected':sequence,'same_protein':sequence},root/'requests')
    request=json.loads((root/'requests'/rows[0]['server_input_file']).read_text())
    # Real Server expands the request version/defaults on return.
    request[0]['version']=3
    request[0]['sequences'][0]['proteinChain']['useStructureTemplate']=True
    cif=gemmi.read_structure(str(pdb(root,sequence))).make_mmcif_document().as_string()
    files={'fold_test_job_request.json':json.dumps(request),'terms_of_use.md':'Original output terms'}
    for index,score in enumerate([.4,.7,.5]):
        files['fold_test_model_'+str(index)+'.cif']=cif
        files['fold_test_summary_confidences_'+str(index)+'.json']=json.dumps(dict(ranking_score=score,ptm=.6))
        files['fold_test_full_data_'+str(index)+'.json']=json.dumps(dict(pae=[[1.]*len(sequence) for _ in sequence]))
    return rows,files


def archive(root,files):
    path=root/'server.zip'
    with zipfile.ZipFile(path,'w') as zipped:
        for name,data in files.items():zipped.writestr(name,data)
    return path


class ServerPreparation(unittest.TestCase):
    def test_server_limits_and_batch_without_sequence_substitution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows=prepare({'short':'MA','long':'A'*5001,'unknown':'MAXA','ready':'MAAA'},root/'requests',limit=1)
            self.assertEqual([r['status'] for r in rows],['unsupported_sequence']*3+['awaiting_prediction'])
            jobs=json.loads((root/'requests/upload.server.json').read_text())
            self.assertEqual(len(jobs),1);self.assertEqual(jobs[0]['modelSeeds'],['42'])
            self.assertEqual(jobs[0]['sequences'][0]['proteinChain']['sequence'],'MAAA')


@unittest.skipUnless(HAS_GEMMI,'Optional gemmi structure parser not installed')
class ServerResult(unittest.TestCase):
    def test_ranked_roundtrip_retains_raw_all_samples_pae_and_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows,files=fixture(root);source=archive(root,files);original=sha256(source)
            catalogue=import_result(root/'requests/requests.json',source,root/'imported','https://alphafoldserver.com/fold/test')
            receipt=json.loads((root/'imported/receipt.json').read_text())
            self.assertEqual(receipt['ranking_score'],.7);self.assertEqual(receipt['returned_samples'],3)
            self.assertEqual(receipt['server_request_version'],3)
            self.assertEqual(receipt['selected_member'],'fold_test_model_1.cif')
            self.assertEqual(sha256(source),original);self.assertEqual(sha256(root/'imported/source.zip'),original)
            dest=root/'stop4';dest.mkdir()
            structures,_=import_structures(catalogue,{'selected':'MAAAA','same_protein':'MAAAA'},dest)
            self.assertEqual(structures['selected']['provider'],'alphafold_server')
            self.assertEqual(structures['selected']['pae_status'],'available_not_interpreted')
            self.assertEqual((root/'imported/terms_of_use.md').read_text(),'Original output terms')
            for path,digest in json.loads((root/'imported/seal.json').read_text()).items():self.assertEqual(sha256(root/'imported'/path),digest)
            with self.assertRaises(ValueError):import_result(root/'requests/requests.json',source,root/'imported')

    def test_changed_prepared_input_is_rejected_before_creating_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows,files=fixture(root);source=archive(root,files)
            (root/'requests'/rows[0]['server_input_file']).write_text('[]')
            with self.assertRaisesRegex(ValueError,'changed'):import_result(root/'requests/requests.json',source,root/'out')
            self.assertFalse((root/'out').exists())

    def test_foreign_sequence_name_or_seed_cannot_be_attached_to_selected_model(self):
        for field in ['sequence','name','seed']:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);rows,files=fixture(root);job=json.loads(files['fold_test_job_request.json'])
                if field=='sequence':job[0]['sequences'][0]['proteinChain']['sequence']='MGGGG'
                if field=='name':job[0]['name']='foreign'
                if field=='seed':job[0]['modelSeeds']=['43']
                files['fold_test_job_request.json']=json.dumps(job)
                with self.assertRaises(ValueError):import_result(root/'requests/requests.json',archive(root,files),root/'out')
                self.assertFalse((root/'out').exists())

    def test_mismatched_coordinate_sequence_has_explicit_failed_attempt(self):
        import gemmi
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows,files=fixture(root)
            files['fold_test_model_1.cif']=gemmi.read_structure(str(pdb(root,'MGGGG'))).make_mmcif_document().as_string()
            with self.assertRaisesRegex(ValueError,'exact requested'):import_result(root/'requests/requests.json',archive(root,files),root/'out')
            self.assertTrue((root/'out/import.failed.json').exists());self.assertFalse((root/'out/receipt.json').exists())

    def test_archive_traversal_is_rejected_without_extracting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows,files=fixture(root);files['../escaped']='no'
            with self.assertRaisesRegex(ValueError,'Unsafe'):import_result(root/'requests/requests.json',archive(root,files),root/'out')
            self.assertFalse((root/'out').exists());self.assertFalse((root/'escaped').exists())

    def test_malformed_pae_fails_instead_of_claiming_ready(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rows,files=fixture(root);files['fold_test_full_data_1.json']=json.dumps(dict(pae=[[1.]]))
            with self.assertRaisesRegex(ValueError,'PAE dimensions'):import_result(root/'requests/requests.json',archive(root,files),root/'out')
            self.assertFalse((root/'out/receipt.json').exists());self.assertTrue((root/'out/import.failed.json').exists())


if __name__=='__main__':unittest.main()
