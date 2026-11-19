import hashlib
from pathlib import Path
import tempfile
import unittest
from genecontext.blast_alignment import alignment_fields
from genecontext.stage1_proteins import parse_protein_hits


class SavedBlastAlignment(unittest.TestCase):
    def fields(self):
        return 'q\ts\t75.000\t4\t1\t3\t2\t5\t1e-8\t30\t3\t5\tA-CD\tABCD'.split('\t')

    def test_saved_trace_and_legacy_keep_same_hsp_id(self):
        fields=self.fields()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'hits.tsv';p.write_text('\t'.join(fields)+'\n')
            h=parse_protein_hits(p,{'q':{'sequence':'ACD'}},{'s':5})[0]
            p.write_text('\t'.join(fields[:12])+'\n')
            old=parse_protein_hits(p,{'q':{'sequence':'ACD'}},{'s':5})[0]
        self.assertEqual(h['hsp_id'],old['hsp_id'])
        self.assertEqual(h['alignment']['qseq'],'A-CD')
        self.assertEqual(h['alignment']['subject_start'],2)
        self.assertNotIn('alignment',old)

    def test_bad_trace_cannot_be_displayed_as_blast(self):
        for index,value in [(12,'A-C'),(13,'ABC'),(12,'X-CD'),(2,'50.000'),(5,'2')]:
            f=self.fields();f[index]=value
            with self.subTest(index=index),self.assertRaises(ValueError):
                alignment_fields(f,'ACD')

    def test_translated_reverse_subject_counts_three_bases_per_residue(self):
        f=self.fields();f[6]='114';f[7]='103'
        a=alignment_fields(f,'ACD',subject_step=3)['alignment']
        self.assertEqual((a['subject_start'],a['subject_end'],a['subject_step']),(114,103,3))


if __name__=='__main__':unittest.main()
