"""Real annotated discovery on unrelated synthetic families, without GFF extents."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from genecontext.project import Project
from genecontext.project_run import run_stage
from genecontext.sequence import sha256

@unittest.skipUnless(all(shutil.which(t) for t in ['blastp','makeblastdb','samtools']),'requires BLAST+ and samtools')
class MinimalStudy(unittest.TestCase):
    def test_protein_and_cds_queries_without_taxon_or_gene_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'study';shutil.copytree(Path(__file__).resolve().parents[1]/'examples/minimal',root)
            inputs={p.name:sha256(p) for p in root.iterdir() if p.is_file()};p=Project(root)
            for target in ['family_a','family_b']:
                result=run_stage(p,target,'stage1');self.assertEqual(result['status'],'awaiting_review',result)
                evidence=json.loads((Path(result['output'])/'evidence.json').read_text())
                self.assertFalse(evidence.get('errors',[]));self.assertEqual(len(evidence['loci']),2)
                self.assertEqual({l['contig_length'] for l in evidence['loci']},{6000})
            self.assertEqual({n:sha256(root/n) for n in inputs},inputs)

if __name__=='__main__':unittest.main()
