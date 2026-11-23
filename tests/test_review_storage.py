"""Report hydration preserves exact datasets and classic script boundaries."""
import unittest
from genecontext.review_storage import finalize


class ReviewStorage(unittest.TestCase):
    def test_idempotent_hydration_keeps_json_and_late_script_order(self):
        dataset = '<script id="interactiveReviewData" type="application/json">{"exact":"\\u003cgene\\u003e"}</script>'
        page = '<html><head></head><body><main></main>'+dataset+'<script>const DATA={};</script></body></html>'
        first = finalize(page)
        self.assertEqual(finalize(first), first)
        late = first.replace('</body>', '<script>window.late=DATA;</script></body>')
        result = finalize(late)
        self.assertIn(dataset, result)
        self.assertEqual(result.count('id="gcfReviewStorageLoader"'), 1)
        self.assertEqual(result.count('// gcf-review-storage-v2'), 1)
        self.assertLess(result.index('const DATA={};'), result.index('window.late=DATA;'))
        self.assertLess(result.index('window.late=DATA;'), result.index('id="gcfReviewStorageLoader"'))
        self.assertIn('type="application/x-gcf-review-script">window.late=DATA;', result)

    def test_noninteractive_reports_are_unchanged(self):
        page = '<html><body>Ordinary report</body></html>'
        self.assertEqual(finalize(page), page)


if __name__ == '__main__':
    unittest.main()
