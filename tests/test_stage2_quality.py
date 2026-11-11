import copy
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.stage2_quality import model_quality, discrepancy_flags, collect_annotation_quality, with_annotation_quality
from genecontext.sequence import sha256


def model(mid, source, blocks, gene='g1', coverage=.8, flags=()):
    return dict(model_id=mid, region_id='r', genome_id='a', gene_group='G',
        original_gene=gene, original_model=mid, source=source, cds_blocks=blocks,
        cds_sha256=mid, target_support=[dict(coverage_subject=coverage)],
        annotation_flags=list(flags))


class ModelRecoveryQuality(unittest.TestCase):
    def split_data(self):
        a=model('a','annotation',[[0,705,0],[1000,1112,0],[2000,2131,2],[3000,3165,0],[4000,4276,0]])
        x=model('x','tiberius',[[0,705,0],[774,783,0]],coverage=0)
        y=model('y','tiberius',[[1779,1783,0],[2000,2131,2],[3000,3165,0],[4000,4276,0]],gene='g2',coverage=.46)
        return dict(models=[a,x,y])

    def test_separate_orfs_jointly_recover_one_annotation_and_remain_separate(self):
        data=self.split_data();before=copy.deepcopy(data);q=model_quality(data)
        self.assertEqual(q['x']['split_models'],['x','y'])
        self.assertTrue(any('split prediction' in w for w in q['y']['warnings']))
        self.assertAlmostEqual(q['x']['annotation_cds_fraction'],705/1389)
        self.assertTrue(any('conflicts' in note for note in q['a']['notes']))
        self.assertTrue(any('Conflicting gene boundaries' in note for note in q['a']['warnings']))
        self.assertEqual(data,before)

    def test_isoforms_or_repeated_predictions_are_not_a_split(self):
        data=self.split_data();data['models'][2]['original_gene']='g1'
        self.assertEqual(model_quality(data)['a']['split_models'],[])
        data=self.split_data()
        for m in data['models'][1:]:m['prediction_model_id']='same_observation'
        self.assertEqual(model_quality(data)['a']['split_models'],[])

    def test_exact_annotation_match_is_not_called_incomplete_from_divergent_query(self):
        data=self.split_data();a=data['models'][0];t=copy.deepcopy(a)
        t.update(model_id='t',source='tiberius',target_support=[dict(coverage_subject=.2)])
        q=model_quality(dict(models=[a,t]))['t']
        self.assertTrue(q['annotation_exact']);self.assertFalse(q['warnings'])

    def test_unannotated_short_query_match_requires_attention(self):
        q=model_quality(dict(models=[model('t','tiberius',[[0,324,0]],coverage=.18)]))['t']
        self.assertTrue(any('below 50%' in w for w in q['warnings']))

    def test_evidence_marker_is_not_a_translation_discrepancy(self):
        self.assertFalse(discrepancy_flags(['exception=annotated by transcript or proteomic data']))
        self.assertTrue(discrepancy_flags(['exception=unclassified translation discrepancy']))

    def test_hash_bound_metadata_retains_transcript_and_cds_flags(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);g=p/'a.gff3';g.write_text(
                'chr\tx\tgene\t1\t30\t.\t+\t.\tID=gene-g\n'
                'chr\tx\tmRNA\t1\t30\t.\t+\t.\tID=rna-r;Parent=gene-g;exception=unclassified transcription discrepancy\n'
                'chr\tx\tCDS\t1\t30\t.\t+\t0\tParent=rna-r;exception=unclassified translation discrepancy\n')
            (p/'a.fna').write_text('>chr\n'+'A'*30+'\n')
            (p/'datasets.tsv').write_text('genome_id\tfasta\tgff3\na\ta.fna\ta.gff3\n')
            m=model('a','annotation',[[0,30,0]]);m.update(original_gene='gene-g',original_model='rna-r')
            data=dict(models=[m],provenance=dict(genomes=[dict(genome_id='a',annotation_sha256=sha256(g))]))
            (p/'evidence.json').write_text(json.dumps(data))
            context=collect_annotation_quality(data,p/'datasets.tsv',sha256(p/'evidence.json'))
            self.assertEqual(len(context['models']['a']),2)
            (p/'annotation_quality.json').write_text(json.dumps(context))
            loaded=with_annotation_quality(p,data)
            self.assertEqual(len(loaded['models'][0]['annotation_flags']),2)
            self.assertFalse(data['models'][0]['annotation_flags'])
            context['source_evidence_sha256']='wrong';(p/'annotation_quality.json').write_text(json.dumps(context))
            with self.assertRaisesRegex(ValueError,'another checkpoint'):with_annotation_quality(p,data)


if __name__=='__main__':unittest.main()
