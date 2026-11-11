import csv
import functools
import hashlib
import http.server
import json
import tempfile
import threading
import unittest
from pathlib import Path
from types import SimpleNamespace

from genecontext.download import checksum_entries, load_sources, run_download, transfer


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.web = self.root / "web"
        self.web.mkdir()
        self.fasta = b">chr1\nATGGCTTTTTAA\n"
        self.gff = b"##gff-version 3\nchr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=tx\n"
        (self.web / "genomic.fna").write_bytes(self.fasta)
        (self.web / "genomic.gff").write_bytes(self.gff)
        (self.web / "md5checksums.txt").write_text("\n".join(
            hashlib.md5(data).hexdigest()+"  ./"+name for name, data in
            [("genomic.fna", self.fasta), ("genomic.gff", self.gff)])+"\n")
        handler = functools.partial(QuietHandler, directory=str(self.web))
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.base = f"http://127.0.0.1:{self.server.server_port}"
        self.manifest = self.root / "sources.tsv"
        self.manifest.write_text("genome_id\taccession\tgenome_url\tgff_url\n"
            f"with_annotation\tGCF_123.1\t{self.base}/genomic.fna\t{self.base}/genomic.gff\n"
            f"no_annotation\tGCA_456.1\t{self.base}/genomic.fna\t\n")

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def args(self, **kwargs):
        return SimpleNamespace(manifest=self.manifest, out=self.root/"downloaded", only=None,
                               include_reference=False, dry_run=False, **kwargs)

    def test_download_checksums_and_generated_manifest(self):
        args = self.args()
        self.assertEqual(run_download(args), 0)
        with open(Path(args.out)/"genomes.tsv") as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1]["gff3"], "")
        for row in rows:
            self.assertEqual((Path(args.out)/row["fasta"]).read_bytes(), self.fasta)
        receipt = json.loads((Path(args.out)/"download_receipt.json").read_text())
        self.assertEqual(receipt["ready"], 2)
        self.assertTrue(all(f["status"] == "downloaded_verified" for r in receipt["genomes"] for f in r["files"]))
        self.assertEqual(run_download(args), 0)
        receipt = json.loads((Path(args.out)/"download_receipt.json").read_text())
        self.assertTrue(all(f["status"] == "reused_verified" for r in receipt["genomes"] for f in r["files"]))

    def test_corrupted_existing_file_preserved_and_next_genome_processed(self):
        args = self.args()
        dest = Path(args.out)/"with_annotation/genome.fna"
        dest.parent.mkdir(parents=True)
        dest.write_bytes(b"changed by user")
        self.assertEqual(run_download(args), 1)
        self.assertEqual(dest.read_bytes(), b"changed by user")
        receipt = json.loads((Path(args.out)/"download_receipt.json").read_text())
        self.assertEqual(receipt["ready"], 1)
        self.assertEqual(receipt["genomes"][0]["status"], "error")
        self.assertEqual(receipt["genomes"][1]["status"], "ready")

    def test_partial_transfer_restart_when_server_ignores_range(self):
        dest = self.root/"genome.fna"
        dest.with_name("genome.fna.part").write_bytes(self.fasta[:7])
        transfer(self.base+"/genomic.fna", dest, hashlib.md5(self.fasta).hexdigest())
        self.assertEqual(dest.read_bytes(), self.fasta)
        self.assertFalse(dest.with_name("genome.fna.part").exists())

    def test_wrong_checksum_never_published_as_ready_file(self):
        dest = self.root/"wrong.fna"
        with self.assertRaises(ValueError):
            transfer(self.base+"/genomic.fna", dest, "0"*32)
        self.assertFalse(dest.exists())
        self.assertTrue(dest.with_name("wrong.fna.part").exists())

    def test_selection_and_invalid_selection(self):
        self.assertEqual(len(load_sources(self.manifest, ["GCF_123.1"])), 1)
        with self.assertRaises(ValueError):
            load_sources(self.manifest, ["wrong_accession"])
        self.assertEqual(checksum_entries("a"*32+"  ./x.gz\n"), {"x.gz":"a"*32})

    def test_ncbi_nested_checksum_basenames_are_distinct(self):
        text = "a"*32+"  ./assembly_structure/Primary_Assembly/chr2acc\n"
        text += "b"*32+"  ./assembly_structure/non-nuclear/chr2acc\n"
        text += "c"*32+"  ./genomic.fna.gz\n"
        parsed = checksum_entries(text)
        self.assertEqual(len(parsed), 3)
        self.assertEqual(parsed["genomic.fna.gz"], "c"*32)

    def test_annotation_only_does_not_fetch_genome_and_keeps_protein_input(self):
        protein=b">NP_1.1 protein\nMAAA\n"
        (self.web/'proteins.faa').write_bytes(protein)
        with (self.web/'md5checksums.txt').open('a') as handle:handle.write(hashlib.md5(protein).hexdigest()+"  ./proteins.faa\n")
        self.manifest.write_text('genome_id\taccession\tspecies\tgff_url\tprotein_url\n'+f"source\tGCF_123.1\tSpecies one\t{self.base}/genomic.gff\t{self.base}/proteins.faa\n")
        args=self.args(annotation_only=True);args.include_reference=True
        self.assertEqual(run_download(args),0)
        with (Path(args.out)/'genomes.tsv').open() as handle:row=next(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(row['fasta'],'');self.assertEqual(row['species'],'Species one');self.assertEqual(row['assembly_accession'],'GCF_123.1')
        self.assertEqual((Path(args.out)/row['protein']).read_bytes(),protein)
        receipt=json.loads((Path(args.out)/'download_receipt.json').read_text());self.assertEqual({f['kind'] for f in receipt['genomes'][0]['files']},{'gff','protein'})

    def test_dry_run_does_not_create_output(self):
        args = self.args()
        args.dry_run = True
        self.assertEqual(run_download(args), 0)
        self.assertFalse(Path(args.out).exists())


if __name__ == "__main__":
    unittest.main()
