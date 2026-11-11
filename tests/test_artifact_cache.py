"""Interrupted and simultaneous builders cannot publish partial artifacts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import threading
import time
import unittest
from genecontext.artifact_cache import cached


class ArtifactCache(unittest.TestCase):
    def test_interruption_does_not_publish_and_retry_succeeds(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            def interrupted(folder):
                (folder/'partial').write_text('unfinished')
                raise KeyboardInterrupt
            with self.assertRaises(KeyboardInterrupt):cached(root,{'input':'same'},interrupted)
            self.assertFalse(list(root.glob('*/cache.json')))
            self.assertFalse(list(root.glob('*.building-*')))
            def complete(folder):(folder/'complete').write_text('finished')
            stored,receipt=cached(root,{'input':'same'},complete)
            self.assertFalse(receipt['reused']);self.assertEqual((stored/'complete').read_text(),'finished')

    def test_concurrent_requests_build_once_and_both_verify_the_result(self):
        with tempfile.TemporaryDirectory() as d:
            barrier=threading.Barrier(2);built=[]
            def build(folder):
                built.append(1);time.sleep(.02);(folder/'result').write_text('complete')
            def request():
                barrier.wait();return cached(Path(d),{'input':'shared'},build)
            with ThreadPoolExecutor(max_workers=2) as pool:
                results=list(pool.map(lambda _:request(),range(2)))
            self.assertEqual(len(built),1)
            self.assertEqual(sorted(receipt['reused'] for _,receipt in results),[False,True])
            self.assertEqual(results[0][0],results[1][0])
