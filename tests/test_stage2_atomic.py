"""Final publication cannot leave a half decision template or half seal."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from genecontext.stage2_resume import atomic_json, pending_table


class AtomicStage2(unittest.TestCase):
    def test_interrupted_decision_write_is_not_visible_and_human_edits_win(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'choices.tsv'
            def interrupted(temporary,*args):
                temporary.write_text('partial');raise KeyboardInterrupt
            with patch('genecontext.sequence.write_tsv',side_effect=interrupted):
                with self.assertRaises(KeyboardInterrupt):pending_table(path,[],['decision'])
            self.assertFalse(path.exists());self.assertEqual(list(Path(d).iterdir()),[])
            pending_table(path,[{'decision':'pending'}],['decision'])
            path.write_text('human choice')
            pending_table(path,[],['decision']);self.assertEqual(path.read_text(),'human choice')

    def test_interrupted_seal_replace_retains_the_previous_complete_file(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'seal.json';atomic_json(path,{'old':'complete'});before=path.read_bytes()
            with patch.object(Path,'replace',side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt):atomic_json(path,{'new':'complete'})
            self.assertEqual(path.read_bytes(),before);self.assertEqual(list(Path(d).iterdir()),[path])
