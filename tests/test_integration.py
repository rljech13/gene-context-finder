"""Real executable integration tests; no search tools are mocked."""
import csv
import gzip
import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from genecontext.sequence import CODE, read_fasta, reverse_complement, translate, write_fasta


TOOLS = ["samtools", "miniprot", "blastn", "blastp", "makeblastdb", "mafft"]
AVAILABLE = all(shutil.which(tool) for tool in TOOLS)


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@unittest.skipUnless(AVAILABLE, "requires real samtools, miniprot, BLAST+ and MAFFT on PATH")
class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        rng = random.Random(71)
        codons = {}
        for codon, aa in CODE.items():
            codons.setdefault(aa, []).append(codon)
        protein = "M" + "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=179))
        self.protein = protein
        self.cds = "".join(rng.choice(codons[aa]) for aa in protein) + "TAA"
        random_dna = lambda n: "".join(rng.choices("ACGT", k=n))
        self.genomic_gene = self.cds[:274] + "GT" + random_dna(96) + "AG" + self.cds[274:]
        prefix, middle = random_dna(250), random_dna(1300)
        minus = reverse_complement(self.cds)
        self.minus_start = len(prefix)+len(self.genomic_gene)+len(middle)
        self.genomic = prefix+self.genomic_gene+middle+minus+random_dna(250)
        self.genome = self.root / "genome with spaces.fna"
        write_fasta(self.genome, [("NC_test.1", self.genomic)])
        self.negative = self.root / "negative.fna.gz"
        with gzip.open(self.negative, "wt") as handle:
            handle.write(">negative\n" + random_dna(4000) + "\n")
        self.query = self.root / "query.faa"
        other = "M" + "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=179))
        write_fasta(self.query, [("target", self.protein), ("negative_query", other)])
        self.annotation = self.root / "genes.gff3"
        start, end = 250, 250+len(self.genomic_gene)
        ms, me = self.minus_start, self.minus_start+len(self.cds)
        self.annotation.write_text(
            "##gff-version 3\n"
            f"NC_test.1\ttest\tgene\t{start+1}\t{end}\t.\t+\t.\tID=gene_plus;Name=TARGET_PLUS\n"
            f"NC_test.1\ttest\tmRNA\t{start+1}\t{end}\t.\t+\t.\tID=tx_plus;Parent=gene_plus\n"
            f"NC_test.1\ttest\tCDS\t{start+1}\t{start+274}\t.\t+\t0\tParent=tx_plus\n"
            f"NC_test.1\ttest\tCDS\t{start+375}\t{end}\t.\t+\t2\tParent=tx_plus\n"
            f"NC_test.1\ttest\tgene\t{ms+1}\t{me}\t.\t-\t.\tID=gene_minus;Name=TARGET_MINUS\n"
            f"NC_test.1\ttest\tmRNA\t{ms+1}\t{me}\t.\t-\t.\tID=tx_minus;Parent=gene_minus\n"
            f"NC_test.1\ttest\tCDS\t{ms+1}\t{me}\t.\t-\t0\tParent=tx_minus\n")

    def tearDown(self):
        self.temp.cleanup()

    def run_pipeline(self, query, kind, rows, name="run", extra=(), expected=0):
        manifest = self.root / f"{name}.tsv"
        manifest.write_text("genome_id\tfasta\tgff3\n" + "\n".join("\t".join(row) for row in rows)+"\n")
        output = self.root / name
        command = [sys.executable, "-m", "genecontext", "run", "--query", str(query),
                   "--query-type", kind, "--genomes", str(manifest), "--out", str(output),
                   "--threads", "1", "--flank", "100", "--context", "1000", *extra]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout+result.stderr+
                         ("\n"+ (output/"manifest.json").read_text() if (output/"manifest.json").exists() else ""))
        return output

    def test_protein_introns_duplicates_both_strands_and_all_genomes(self):
        output = self.run_pipeline(self.query, "protein", [
            ("positive", str(self.genome), ""), ("negative", str(self.negative), "")])
        summary = table(output/"genome_query_summary.tsv")
        self.assertEqual(len(summary), 4)
        self.assertEqual([r["status"] for r in summary if r["genome_id"] == "negative"], ["no_hit", "no_hit"])
        models = [r for r in table(output/"models.tsv") if r["query"] == "q0001"]
        complete = [r for r in models if r["status"] == "complete_predicted"]
        self.assertEqual(len(complete), 2, models)
        self.assertEqual({r["strand"] for r in complete}, {"+", "-"})
        self.assertTrue(any(len(json.loads(r["cds_exons"])) == 2 for r in complete))
        extracted = read_fasta(output/"sequences/q0001/cds.fna")
        for row in complete:
            self.assertEqual(extracted[row["model_id"]], self.cds)
        alignment = read_fasta(output/"alignments/q0001/cds.codon.aligned.fna")
        self.assertEqual(len(alignment), 2)
        for aligned in alignment.values():
            self.assertEqual(translate(aligned.replace("-", "")), self.protein)
        self.assertFalse(Path(str(self.genome)+".fai").exists())
        self.assertEqual(len(list((output/"genomes/positive").glob("*.context.svg"))), len(models))
        original_manifest = (output/"manifest.json").read_bytes()
        realigned = self.root/"realigned"
        result = subprocess.run([sys.executable, "-m", "genecontext", "align", "--results", str(output),
            "--out", str(realigned), "--threads", "1", "--dna-alignment", "fast"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((output/"manifest.json").read_bytes(), original_manifest)
        new_manifest = json.loads((realigned/"manifest.json").read_text())
        self.assertEqual(new_manifest["status"], "completed")
        self.assertEqual(new_manifest["source_run_status"], "completed")
        self.assertTrue((realigned/"alignments/q0001/cds.codon.aligned.fna").exists())
        self.assertTrue(any("--maxiterate" in c["argv"] for c in new_manifest["commands"]))

    def test_genomic_dna_recovers_annotation_and_noncoding_stays_dna(self):
        query = self.root / "genomic_query.fna"
        write_fasta(query, [("genomic_gene", self.genomic_gene)])
        output = self.run_pipeline(query, "dna", [("annotated", str(self.genome), str(self.annotation))])
        models = table(output/"models.tsv")
        self.assertEqual(len(models), 2)
        self.assertEqual({r["source_model"] for r in models}, {"tx_plus", "tx_minus"})
        self.assertTrue(all(r["status"] == "complete_predicted" for r in models))
        self.assertEqual(set(read_fasta(output/"sequences/q0001/cds.fna").values()), {self.cds})
        self.assertGreaterEqual(len(table(output/"neighbors.tsv")), 2)
        output2 = self.run_pipeline(query, "dna", [("unannotated", str(self.genome), "")], name="raw_dna")
        self.assertEqual(table(output2/"genome_query_summary.tsv")[0]["status"], "dna_hits_no_cds_model")
        self.assertTrue(table(output2/"dna_regions.tsv"))
        self.assertEqual(table(output2/"models.tsv"), [])
        self.assertTrue(read_fasta(output2/"alignments/q0001/dna_regions.aligned.fa"))
        self.assertTrue(read_fasta(output2/"sequences/q0001/dna_region_upstream.fna"))

    def test_cds_route_and_annotation_deduplication(self):
        query = self.root / "coding_query.fna"
        write_fasta(query, [("coding", self.cds)])
        output = self.run_pipeline(query, "cds", [("annotated", str(self.genome), str(self.annotation))])
        models = table(output/"models.tsv")
        self.assertEqual(len(models), 2, models)
        self.assertTrue(all(r["source"] == "annotation" for r in models))
        codon = read_fasta(output/"alignments/q0001/cds.codon.aligned.fna")
        self.assertEqual(len(codon), 3)  # Two genomic loci and the known query CDS.

    def test_broken_genome_is_reported_and_later_genome_processed(self):
        output = self.run_pipeline(self.query, "protein", [
            ("broken", str(self.root/"missing.fna"), ""),
            ("later", str(self.negative), "")], extra=["--skip-alignments"], expected=1)
        summary = table(output/"genome_query_summary.tsv")
        self.assertEqual(len(summary), 4)
        self.assertEqual([r["status"] for r in summary], ["error", "error", "no_hit", "no_hit"])
        self.assertEqual(json.loads((output/"manifest.json").read_text())["status"], "completed_with_errors")

    def test_exact_assembly_annotation_mismatch_is_error(self):
        bad = self.root / "mismatch.gff3"
        bad.write_text(self.annotation.read_text().replace("NC_test.1", "another_assembly_chr"))
        output = self.run_pipeline(self.query, "protein", [("mismatch", str(self.genome), str(bad))],
                                   extra=["--skip-alignments"], expected=1)
        summary = table(output/"genome_query_summary.tsv")
        self.assertTrue(all(r["status"] == "error" for r in summary))
        self.assertIn("mismatch", summary[0]["error"].lower())

    def test_long_pipe_contig_ids_survive_nucleotide_search(self):
        long_id = "custom_" + "X"*90
        piped = "ref|NC_000001.1|"
        genome = self.root / "long_headers.fna"
        write_fasta(genome, [(long_id, self.genomic_gene), (piped, self.genomic_gene)])
        query = self.root / "genomic_query.fna"
        write_fasta(query, [("genomic", self.genomic_gene)])
        output = self.run_pipeline(query, "dna", [("long_headers", str(genome), "")],
                                   extra=["--skip-alignments"])
        self.assertEqual({r["contig"] for r in table(output/"dna_regions.tsv")}, {long_id, piped})


if __name__ == "__main__":
    unittest.main()
