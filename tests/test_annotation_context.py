import gzip
from pathlib import Path
import tempfile
import unittest

from genecontext.annotation_context import read_features, collect
from genecontext.sequence import sha256


class FullGeneContext(unittest.TestCase):
    def test_noncoding_pseudogene_and_boundary_overlap_are_retained(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'genes.gff.gz'
            with gzip.open(p,'wt') as f:
                f.write('chr\tx\tgene\t11\t30\t.\t+\t.\tID=coding;Name=PROTEIN;gene_biotype=protein_coding\n'
                        'chr\tx\tgene\t41\t70\t.\t-\t.\tID=rna;Name=RNA;gene_biotype=lncRNA\n'
                        'chr\tx\texon\t41\t70\t.\t-\t.\tParent=rna\n'
                        'chr\tx\tpseudogene\t91\t110\t.\t+\t.\tID=pseudo;Name=PSEUDO\n'
                        'chr\tx\tgene\t111\t120\t.\t+\t.\tID=outside\n')
            region=dict(region_id='r',genome_id='g',contig='chr',window_start=10,window_end=110)
            rows=read_features(p,[region]);self.assertEqual([r['original_gene'] for r in rows],['coding','rna','pseudo'])
            self.assertEqual((rows[1]['start'],rows[1]['end'],rows[1]['biotype']),(40,70,'lncRNA'))
            other=dict(region,region_id='second');copies=read_features(p,[region,other])
            self.assertEqual(copies[0]['feature_id'],copies[1]['feature_id'])

    def test_existing_checkpoint_requires_the_exact_original_gff(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);gff=root/'g.gff';gff.write_text('chr\tx\tgene\t1\t30\t.\t+\t.\tID=rna;gene_biotype=lncRNA\n')
            manifest=root/'genomes.tsv';manifest.write_text('genome_id\tfasta\tgff3\ng\tgenome.fa\tg.gff\n')
            data=dict(regions=[dict(region_id='r',genome_id='g',contig='chr',window_start=0,window_end=100)],
                      provenance=dict(genomes=[dict(genome_id='g',annotation_sha256=sha256(gff))]))
            context=collect(data,manifest);self.assertEqual(len(context['features']),1)
            self.assertEqual(context['sources'][0]['sha256'],sha256(gff))
            gff.write_text(gff.read_text()+'#changed release\n')
            with self.assertRaisesRegex(ValueError,'differs from the original'):collect(data,manifest)
