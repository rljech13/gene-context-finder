import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from genecontext.alignment import backtranslate, pairwise_stats
from genecontext.pipeline import extract_model, load_genomes, load_queries, merge_dna_regions
from genecontext.sequence import (Intervals, Model, coding_sequence, read_fasta, read_gff,
                                  reverse_complement, translate, union_length)


class MemoryGenome:
    def __init__(self, sequence):
        self.sequence = sequence

    def length(self, contig):
        return len(self.sequence)

    def fetch(self, contig, start, end, strand="+"):
        sequence = self.sequence[start:end]
        return reverse_complement(sequence) if strand == "-" else sequence


class SequenceTests(unittest.TestCase):
    def test_code_and_reverse_complement(self):
        self.assertEqual(translate("ATGGCTTTCTAA"), "MAF*")
        self.assertEqual(reverse_complement("ACGTNRY"), "RYNACGT")
        self.assertEqual(union_length([(1, 5), (2, 3), (4, 8), (10, 12)]), 9)

    def test_phase_one_intron_keeps_split_codon(self):
        # Four bases before intron: the two remaining bases complete codon GCT.
        genome = MemoryGenome("ATGG" + "GTCCAG" + "CTTTCTAA")
        model = Model("tx", "chr", "+", [(0, 4, 0), (10, 18, 2)])
        cds, protein, status, issues = coding_sequence(model, genome)
        self.assertEqual((cds, protein, status, issues),
                         ("ATGGCTTTCTAA", "MAF", "complete_predicted", []))

    def test_minus_strand_and_separate_stop(self):
        genomic = reverse_complement("ATGGCTTTCTAA")
        genome = MemoryGenome("AAA" + genomic + "GGG")
        model = Model("tx", "chr", "-", [(6, 15, 0)], stops=[(3, 6)])
        self.assertEqual(coding_sequence(model, genome),
                         ("ATGGCTTTCTAA", "MAF", "complete_predicted", []))

    def test_phases_and_disruption_are_not_silently_repaired(self):
        model = Model("bad", "chr", "+", [(0, 4, 0), (4, 12, 1)])
        self.assertIn("inconsistent_phase", coding_sequence(model, MemoryGenome("ATGGCTTTCTAA"))[3])
        model = Model("bad", "chr", "+", [(0, 12, 0)])
        self.assertEqual(coding_sequence(model, MemoryGenome("ATGTAATTTTAA"))[2], "disrupted")
        self.assertEqual(coding_sequence(Model("p", "chr", "+", [(0, 9, 0)]),
                                         MemoryGenome("ATGGCTTTT"))[2], "partial")

    def test_stop_overlapping_cds_not_duplicated(self):
        model = Model("tx", "chr", "+", [(0, 12, 0)], stops=[(9, 12)])
        self.assertEqual(coding_sequence(model, MemoryGenome("ATGGCTTTCTAA"))[0], "ATGGCTTTCTAA")

    def test_flanks_are_oriented_and_clipped(self):
        genome = MemoryGenome("AAA" + reverse_complement("ATGGCTTAA") + "CCCCC")
        model = Model("tx", "chr", "-", [(3, 12, 0)])
        row, seqs = extract_model(SimpleNamespace(flank=5), genome, model, "m", "q", "g",
                                 dict(coverage=1, identity=1, evidence="fixture"))
        self.assertEqual(seqs["upstream"], "GGGGG")
        self.assertEqual(seqs["downstream"], "TTT")
        self.assertTrue(row["upstream_complete"])
        self.assertFalse(row["downstream_complete"])
        self.assertEqual(row["downstream_start"], 0)

    def test_codon_alignment_roundtrip_and_mismatch(self):
        self.assertEqual(backtranslate({"a": "MA-F", "b": "M-GF"},
                                       {"a": "ATGGCTTTT", "b": "ATGGGTTTT"}),
                         {"a": "ATGGCT---TTT", "b": "ATG---GGTTTT"})
        with self.assertRaises(ValueError):
            backtranslate({"a": "MK"}, {"a": "ATGGCT"})

    def test_identity_coverage_excludes_ambiguous_positions(self):
        row = pairwise_stats({"a": "ACTN--", "b": "AC-GTT"}, set("ACGT"))[0]
        self.assertEqual(row["identity"], 1)
        self.assertEqual(row["paired_positions"], 2)
        self.assertEqual(row["coverage_a"], 2/3)
        self.assertEqual(row["coverage_b"], 2/5)

    def test_nested_interval_lookup(self):
        records = [(0, 10000), (5, 10), (20, 30)]
        index = Intervals(records, lambda r: "chr", lambda r: r[0], lambda r: r[1])
        self.assertEqual(index.overlap("chr", 500, 510), [(0, 10000)])


class FileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_shared_parents_and_gz_gff(self):
        path = self.root / "genes.gff3.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("##gff-version 3\nchr\ts\tCDS\t1\t12\t.\t+\t0\tParent=tx1,tx2\n")
        models, _ = read_gff(path)
        self.assertEqual([m.id for m in models], ["tx1", "tx2"])
        self.assertEqual(models[0].cds, [(0, 12, 0)])

    def test_duplicate_query_headers_rejected(self):
        path = self.root / "q.fa"
        path.write_text(">same\nACGT\n>same\nACGT\n")
        with self.assertRaises(ValueError):
            read_fasta(path)

    def test_manifest_paths_relative_and_all_genomes_kept(self):
        path = self.root / "genomes.tsv"
        path.write_text("genome_id\tfasta\tgff3\na\tx.fna\tx.gff3\nb\ty.fa.gz\t\n")
        rows = load_genomes(SimpleNamespace(genomes=path))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["fasta"], str((self.root / "x.fna").resolve()))
        self.assertEqual(rows[1]["gff3"], "")

    def test_directory_recursive_gzip_and_duplicate_stems(self):
        (self.root / "nested").mkdir()
        (self.root / "a.fa").touch()
        (self.root / "nested" / "b.fna.gz").touch()
        args = SimpleNamespace(genomes=None, genome_dir=self.root)
        self.assertEqual(len(load_genomes(args)), 2)
        (self.root / "nested" / "a.fasta").touch()
        with self.assertRaises(ValueError):
            load_genomes(args)

    def test_dna_not_forced_to_code(self):
        path = self.root / "query.fa"
        path.write_text(">region\n" + "ACGT" * 30 + "\n")
        queries = load_queries(SimpleNamespace(query=path, query_type="dna"))
        self.assertFalse(queries[0]["protein"])
        path.write_text(">cds\nATG" + "GCT"*35 + "TAA\n")
        queries = load_queries(SimpleNamespace(query=path, query_type="dna"))
        self.assertTrue(queries[0]["inferred_cds"])

    def test_dna_region_merging_respects_query_strand_and_contig(self):
        hsps = [dict(query="q", contig="chr", strand="+", start=10, end=20),
                dict(query="q", contig="chr", strand="+", start=25, end=35),
                dict(query="q", contig="chr", strand="-", start=15, end=30),
                dict(query="other", contig="chr", strand="+", start=10, end=20)]
        regions = merge_dna_regions(hsps, 10)
        self.assertEqual(len(regions), 3)
        self.assertEqual(max(r["hsp_count"] for r in regions), 2)


if __name__ == "__main__":
    unittest.main()
