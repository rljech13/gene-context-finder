import functools
import gzip
import hashlib
import http.server
import json
from pathlib import Path
import tempfile
import threading
import unittest

from genecontext.context_compare import prune_reference_catalog
from genecontext.context_reference import add_references, swissprot


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


class ContextReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.web = self.root / 'web'; self.web.mkdir()
        self.payload = gzip.compress(b'>sp|TEST|GENE Protein OS=Test GN=GENE\n' + b'MACDEFGHIKL'*8 + b'\n')
        (self.web / 'uniprot_sprot.fasta.gz').write_bytes(self.payload)
        self.expected = hashlib.md5(self.payload).hexdigest()
        self.metalink()
        handler = functools.partial(QuietHandler, directory=str(self.web))
        self.server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True); self.thread.start()
        self.base = f'http://127.0.0.1:{self.server.server_port}/'

    def metalink(self):
        (self.web / 'RELEASE.metalink').write_text(f'''<metalink xmlns="http://www.metalinker.org/" version="3.0"><version>test_release</version>
        <files><file name="uniprot_sprot.fasta.gz"><size>{len(self.payload)}</size><verification>
        <hash type="md5">{self.expected}</hash></verification></file></files></metalink>''')

    def tearDown(self):
        self.server.shutdown(); self.server.server_close(); self.thread.join(); self.temp.cleanup()

    def test_real_http_download_official_checksum_and_offline_reuse(self):
        path, info = swissprot(self.root / 'cache', self.base)
        self.assertEqual(path.read_bytes(), self.payload)
        self.assertEqual(info['release'], 'test_release')
        self.assertEqual(info['sha256'], hashlib.sha256(self.payload).hexdigest())
        # Metadata endpoint disappears: reuse must not contact it again.
        (self.web / 'RELEASE.metalink').unlink()
        self.assertEqual(swissprot(self.root / 'cache', self.base)[1]['cache_status'], 'reused_verified')
        path.write_bytes(b'corrupted')
        with self.assertRaisesRegex(ValueError, 'fails verification'):
            swissprot(self.root / 'cache', self.base)
        self.assertEqual(path.read_bytes(), b'corrupted')
        path.write_bytes(self.payload)
        (path.parent / 'receipt.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'Invalid Swiss-Prot cache receipt'):
            swissprot(self.root / 'cache', self.base)

    def test_wrong_upstream_checksum_never_becomes_ready(self):
        self.expected = '0'*32; self.metalink()
        with self.assertRaisesRegex(ValueError, 'MD5 mismatch'):
            swissprot(self.root / 'cache', self.base)
        self.assertFalse((self.root / 'cache/swissprot/receipt.json').exists())
        self.assertFalse((self.root / 'cache/swissprot/uniprot_sprot.fasta.gz').exists())
        self.assertFalse((self.root / 'cache/swissprot/download.lock').exists())

    def test_reference_names_skips_and_local_catalog_pruning(self):
        path = self.root / 'proteins.faa'
        path.write_text('>sp|ID|ENTRY Full description OS=Test GN=SYMBOL PE=1\n' + 'M'*60 + '\n>short\nMA\n')
        proteins, metadata = {}, {}
        stats = add_references(path, proteins, metadata, 'swissprot_reference', skip_invalid=True)
        self.assertEqual(stats, dict(usable_proteins=1, skipped_short_or_nonstandard=1))
        rid = next(iter(proteins))
        self.assertEqual(metadata[rid]['name'], 'SYMBOL')
        self.assertEqual(metadata[rid]['gene_id'], 'sp|ID|ENTRY')
        proteins['unused'] = 'M'*60; metadata['unused'] = dict(genome_id='external_reference')
        proteins['prediction'] = 'M'*60
        metadata['prediction'] = dict(genome_id='genome', reference_protein=rid)
        self.assertEqual(prune_reference_catalog(proteins, metadata), 1)
        self.assertEqual(set(proteins), {rid, 'prediction'})
        with self.assertRaisesRegex(ValueError, 'Invalid context'):
            add_references(path, {}, {}, 'supplied_reference')


if __name__ == '__main__':
    unittest.main()
