import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.sequence import reverse_complement,sha256
from genecontext.stage2_reconstruct import transcripts,project_sam,load
from genecontext.stage2_quality import model_quality

class TranscriptReconstruction(unittest.TestCase):
    def fixture(self,reverse=False):
        seq='GGATGAAATAACC';ref=dict(accession='RNA.1',sequence=seq,cds_start=2,cds_end=11)
        target=reverse_complement(seq) if reverse else seq
        dna='N'*5+target[:6]+'GTCCCCCCAG'+target[6:]+'N'*5
        sam=f'RNA.1\t{16 if reverse else 0}\tr\t6\t60\t6M10N7M\t*\t0\t0\t{target}\t*\tNM:i:0\n'
        return ref,dna,sam

    def test_exact_splice_projection_preserves_split_codons_and_utr(self):
        ref,dna,sam=self.fixture();r=project_sam(sam,ref,'r',dna,100)
        self.assertEqual(r['cds'],'ATGAAATAA');self.assertEqual(r['protein'],'MK')
        self.assertEqual(r['cds_blocks'],[[107,111,0],[121,126,2]])
        self.assertEqual(r['transcript_length'],13)

    def test_reverse_strand_maps_original_transcript_cds_and_phases(self):
        ref,dna,sam=self.fixture(True);r=project_sam(sam,ref,'r',dna,100)
        self.assertEqual(r['strand'],'-');self.assertEqual(r['cds'],'ATGAAATAA')
        self.assertEqual(r['cds_blocks'],[[121,126,0],[107,111,1]])

    def test_bad_mapping_is_not_admitted_as_repaired_gene(self):
        ref,dna,sam=self.fixture()
        for bad in [sam.replace('\t60\t','\t0\t'),sam.replace('NM:i:0','NM:i:1'),sam.replace('6M10N7M','1S5M10N7M')]:
            with self.assertRaises(ValueError):project_sam(bad,ref,'r',dna,100)
        with self.assertRaises(ValueError):project_sam(sam,ref,'r',dna.replace('ATG','ACG'),100)

    def test_model_support_is_bound_to_same_checkpoint_and_same_cds(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);ref,dna,sam=self.fixture();rec=project_sam(sam,ref,'r',dna,100)
            mid='m';h=lambda s:hashlib.sha256(s.encode()).hexdigest()
            model=dict(model_id=mid,source='annotation',cds_blocks=copy.deepcopy(rec['cds_blocks']),strand='+',cds_sha256=h(rec['cds']),protein_sha256=h(rec['protein']))
            (p/'evidence.json').write_text(json.dumps(dict(models=[model])))
            rec.update(model_id=mid,cds_sha256=model['cds_sha256'],protein_sha256=model['protein_sha256'],status='existing_full_model_reconstructed',rna_all_introns_reported=True,resolves_boundary_conflict=True)
            payload=dict(schema='stage2_reconstructions_v1',source_evidence_sha256=sha256(p/'evidence.json'),models=[rec]);(p/'reconstructions.json').write_text(json.dumps(payload))
            data=dict(models=[model]);result=load(p,data)
            self.assertTrue(result['models'][0]['reconstruction']['resolves_boundary_conflict'])
            self.assertNotIn('reconstruction',model)
            payload['models'][0]['cds_blocks'][0][0]+=1;(p/'reconstructions.json').write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'differs'):load(p,data)
            payload['source_evidence_sha256']='other';(p/'reconstructions.json').write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'another checkpoint'):load(p,data)

    def test_predicted_rna_alone_does_not_resolve_connectivity_conflict(self):
        from tests.test_stage2_quality import ModelRecoveryQuality
        data=ModelRecoveryQuality().split_data();a=data['models'][0]
        a['reconstruction']=dict(accession='RNA.1',resolves_boundary_conflict=False)
        self.assertTrue(model_quality(data)['a']['warnings'])
        a['reconstruction']['resolves_boundary_conflict']=True
        q=model_quality(data)
        self.assertEqual(q['a']['warnings'],[])
        self.assertTrue(q['x']['warnings']);self.assertTrue(q['y']['warnings'])

    def test_genbank_parser_retains_reported_rna_evidence_without_claiming_own_experiment(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'rna.xml';p.write_text('<GBSet><GBSeq><GBSeq_accession-version>RNA.1</GBSeq_accession-version><GBSeq_definition>PREDICTED mRNA</GBSeq_definition><GBSeq_sequence>GGATGAAATAACC</GBSeq_sequence><GBSeq_feature-table><GBFeature><GBFeature_key>gene</GBFeature_key><GBFeature_quals><GBQualifier><GBQualifier_name>note</GBQualifier_name><GBQualifier_value>RNAseq alignments, support for all annotated introns</GBQualifier_value></GBQualifier></GBFeature_quals></GBFeature><GBFeature><GBFeature_key>CDS</GBFeature_key><GBFeature_location>3..11</GBFeature_location></GBFeature></GBSeq_feature-table></GBSeq></GBSet>')
            r=transcripts(p)['RNA.1'];self.assertTrue(r['reference_derived']);self.assertTrue(r['rna_all_introns_reported']);self.assertEqual((r['cds_start'],r['cds_end']),(2,11))
