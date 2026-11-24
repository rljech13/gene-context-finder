"""Exact-input confirmation gates automatic core results, never client numbers."""
from copy import deepcopy
import hashlib
from http.client import HTTPConnection
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext import tree_placement as core
from genecontext.tree_review import write_bundle, IQTREE_DOWNLOADS
from genecontext.tree_review_service import (validate_panel, prepare_overlay, default_review,
    create_server, ReviewState, orient_tree, rebuild, alignment_input, REVIEW_SCHEMA)


def digest(x): return hashlib.sha256(x.encode()).hexdigest()


def panels():
    rows = {n: 'MADCEFGHIKLMNPQRSTVWY' * 2 for n in 'abcde'}
    models = [{'model_id': 'model_' + n.upper(), 'tree_id': n, 'taxon': n.upper(),
        'primary': True, 'protein_sha256': digest(rows[n]), 'label': 'Organism ' + n} for n in 'abcde']
    leaves = [{'tree_id': n, 'sequence_ids': ['model_' + n.upper()], 'labels': ['Organism ' + n],
        'protein_sha256': digest(rows[n]), 'included': n != 'e'} for n in rows]
    return {'test': {'seal': 'sealed-source', 'label': 'Any family', 'groups': {'family': {
        'label': 'Family', 'rows': rows, 'models': models, 'leaves': leaves,
        'trees': {'saved': '((a:1,b:1)80:1,(c:1,d:1)80:1);', 'other': '((a:1,c:1)100:1,(b:1,d:1)100:1);'},
        'species_tree': '((A,B),(C,D),E);', 'inputSha256': 'stable-input'}}}}


class ReviewContract(unittest.TestCase):
    def setUp(self):
        self.manifest = {'schema': 'gcf_tree_review_inputs_v2', 'panels': panels()}
        self.reviews = default_review(self.manifest)
        self.item = self.reviews['panels']['test']
        self.request = self.item['groups']['family']

    def approved(self):
        preview = validate_panel(self.manifest, self.item)
        self.request.update(approved=True, confirmed_signature=preview['groups']['family']['signature'])
        return validate_panel(self.manifest, self.item)

    def test_pending_is_proposed_only_then_exact_confirmation_enables(self):
        result = validate_panel(self.manifest, self.item)
        for n in 'ABCD':
            self.assertEqual(result['models']['model_' + n]['proposed_points'], 12)
            self.assertIsNone(result['models']['model_' + n]['points'])
        approved = self.approved()
        self.assertTrue(approved['groups']['family']['approved'])
        self.assertEqual(approved['models']['model_A']['points'], 12)

    def test_default_membership_uses_actual_tree_tips_not_rows(self):
        self.assertEqual(self.request['selected_ids'], list('abcd'))
        self.assertIn('e', self.manifest['panels']['test']['groups']['family']['rows'])

    def test_stage3_mask_and_outgroup_are_preserved_but_never_preapproved(self):
        source = self.manifest['panels']['test']['groups']['family']
        source.update(mask=[0, 1, 2, 3], outgroup_ids=['a'], root_basis='Stage3 choice')
        item = default_review(self.manifest)['panels']['test']
        group = validate_panel(self.manifest, item)['groups']['family']
        self.assertEqual(group['outgroup_ids'], ['a']); self.assertEqual(group['mask'], [0, 1, 2, 3])
        self.assertFalse(group['approved']); self.assertFalse(group['requires_rebuild'])

    def test_omitted_focal_model_is_na_even_after_review(self):
        result = self.approved()['models']['model_E']
        self.assertIsNone(result['points']); self.assertIsNone(result['proposed_points'])
        self.assertEqual(result['reason'], 'model_not_in_tree')

    def test_outgroup_change_clears_approval_but_preserves_proposed_points(self):
        first = self.approved()
        self.request.update(outgroup_ids=['a'], root_basis='Declared outgroup')
        changed = validate_panel(self.manifest, self.item)
        self.assertFalse(changed['groups']['family']['approved'])
        self.assertNotEqual(first['groups']['family']['signature'], changed['groups']['family']['signature'])
        for mid in first['models']:
            self.assertEqual(first['models'][mid]['proposed_points'], changed['models'][mid]['proposed_points'])
        self.assertNotEqual(changed['groups']['family']['newick'], changed['groups']['family']['display_newick'])

    def test_explicitly_unrooted_tree_can_be_approved(self):
        result = self.approved()['groups']['family']
        self.assertEqual(result['outgroup_ids'], [])
        self.assertEqual(result['newick'], result['display_newick'])

    def test_tree_switch_invalidates_and_recomputes(self):
        self.approved(); self.request['source'] = {'id': 'other'}
        result = validate_panel(self.manifest, self.item)
        self.assertFalse(result['groups']['family']['approved'])
        self.assertEqual(result['models']['model_A']['proposed_points'], 0)

    def test_mask_and_membership_and_build_method_require_rebuild(self):
        for change in [{'mask': [0, 1, 2, 3]}, {'selected_ids': list('abc')}, {'tree_method': {'model': 'LG'}}]:
            with self.subTest(change=change):
                self.setUp(); self.approved(); self.request.update(change)
                result = validate_panel(self.manifest, self.item)['groups']['family']
                self.assertFalse(result['approved']); self.assertTrue(result['requires_rebuild'])

    def test_reference_and_core_method_options_invalidate_confirmation(self):
        self.approved()
        self.manifest['panels']['test']['groups']['family']['species_tree'] = '((A,C),(B,D),E);'
        result = validate_panel(self.manifest, self.item)
        self.assertFalse(result['groups']['family']['approved'])
        self.assertEqual(result['models']['model_A']['proposed_points'], 0)
        self.setUp(); self.approved(); self.request['method'] = {'support_field': 0}
        self.assertFalse(validate_panel(self.manifest, self.item)['groups']['family']['approved'])

    def test_model_sequence_and_alignment_pool_invalidate_confirmation(self):
        self.approved(); group = self.manifest['panels']['test']['groups']['family']
        group['models'][0]['protein_sha256'] = digest('different protein')
        self.assertFalse(validate_panel(self.manifest, self.item)['groups']['family']['approved'])
        self.setUp(); self.approved(); rows = self.manifest['panels']['test']['groups']['family']['rows']; rows['a'] = 'A' * len(rows['a'])
        self.assertFalse(validate_panel(self.manifest, self.item)['groups']['family']['approved'])

    def test_core_source_bytes_bound_to_confirmation(self):
        self.approved()
        with patch('genecontext.tree_review_service.core.METHOD', 'new_core_method'):
            self.assertFalse(validate_panel(self.manifest, self.item)['groups']['family']['approved'])

    def test_forged_numeric_points_split_and_decisions_are_ignored(self):
        self.item['models'] = {'model_A': {'points': 15}}
        self.request.update(points=15, decisions={'model_A': {'points': 15}}, split={'support': 100},
            placement={'models': {'model_A': {'points': 15}}})
        result = self.approved()
        self.assertEqual(result['models']['model_A']['points'], 12)
        self.assertNotIn('decisions', result['request']['groups']['family'])

    def test_checkbox_without_matching_preview_is_not_approval(self):
        for value in [None, 'invented-signature']:
            self.request.update(approved=True, confirmed_signature=value)
            self.assertFalse(validate_panel(self.manifest, self.item)['groups']['family']['approved'])

    def test_request_roundtrip_revalidates_approved_overlay(self):
        approved = self.approved()
        overlay = prepare_overlay({'schema': REVIEW_SCHEMA, 'panels': {'test': approved['request']}}, self.manifest)
        self.assertEqual(overlay['schema'], 'gcf_tree_review_overlay_v2')
        self.assertEqual(overlay['panels']['test']['models']['model_A']['points'], 12)
        self.assertEqual(prepare_overlay(overlay['request'], self.manifest)['panels']['test']['models']['model_A']['points'], 12)

    def test_wrong_panel_input_and_old_manual_schema_rejected(self):
        for mutation in [{'seal': 'other-seal'}, {'panel': 'unknown'}]:
            with self.assertRaises(ValueError): validate_panel(self.manifest, dict(self.item, **mutation))
        self.request['inputSha256'] = 'wrong'
        with self.assertRaises(ValueError): validate_panel(self.manifest, self.item)
        with self.assertRaises(ValueError): prepare_overlay({'schema': 'gcf_tree_review_v1', 'panels': {}}, self.manifest)

    def test_cannot_replace_saved_tree_under_its_id_or_checksum(self):
        self.request['source']['sha256'] = 'fake'
        with self.assertRaises(ValueError): validate_panel(self.manifest, self.item)
        self.request['source'] = {'id': 'saved', 'newick': '((a,c)100,(b,d)100);'}
        self.assertEqual(validate_panel(self.manifest, self.item)['models']['model_A']['proposed_points'], 12)

    def test_outgroup_must_be_monophyletic_and_have_basis(self):
        self.request.update(outgroup_ids=['a', 'c'], root_basis='Proposed')
        with self.assertRaisesRegex(ValueError, 'not monophyletic'): validate_panel(self.manifest, self.item)
        self.request.update(outgroup_ids=['a'], root_basis='')
        with self.assertRaisesRegex(ValueError, 'basis'): validate_panel(self.manifest, self.item)

    def test_identical_aliases_cannot_be_divided(self):
        self.manifest['panels']['test']['groups']['family']['leaves'][0]['sequence_ids'].append('reference_A')
        self.request.update(outgroup_ids=['model_A'], root_basis='Proposed')
        with self.assertRaisesRegex(ValueError, 'identical aliases'): validate_panel(self.manifest, self.item)
        self.request['outgroup_ids'] = ['model_A', 'reference_A']
        self.assertIn('display_newick', validate_panel(self.manifest, self.item)['groups']['family'])

    def test_root_orientation_preserves_core_quartets(self):
        original = self.manifest['panels']['test']['groups']['family']['trees']['saved']
        models = self.manifest['panels']['test']['groups']['family']['models']
        reference = self.manifest['panels']['test']['groups']['family']['species_tree']
        first = core.assess(original, models, species_tree=reference)
        for outgroup in [['a'], ['a', 'b'], ['c', 'd']]:
            rerooted = orient_tree(original, outgroup)
            other = core.assess(rerooted, models, species_tree=reference)
            self.assertEqual({k: v['points'] for k, v in first['models'].items()}, {k: v['points'] for k, v in other['models'].items()})

    def test_rebuilt_source_bound_to_exact_alignment_and_rows(self):
        group = self.manifest['panels']['test']['groups']['family']
        selected, mask, columns, fasta = alignment_input(group, list('abcd'))
        text = group['trees']['other']
        self.request.update(source={'id': 'rebuilt:test', 'sha256': digest(text), 'newick': text,
            'provenance': {'kind': 'rebuilt', 'alignment_sha256': digest(fasta), 'selected_ids': selected, 'mask': mask, 'columns': columns}})
        self.assertEqual(self.approved()['models']['model_A']['points'], 0)
        self.request['source']['provenance']['alignment_sha256'] = 'forged'
        self.assertTrue(validate_panel(self.manifest, self.item)['groups']['family']['requires_rebuild'])


class PortableBundle(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.manifest = write_bundle(self.root, panels())
        self.source = Path(__file__).resolve().parents[1] / 'genecontext'

    def tearDown(self): self.tmp.cleanup()

    def test_exact_ui_core_and_server_bytes_ship(self):
        for source, target in [('tree_placement.py', 'tree_placement.py'), ('tree_review_service.py', 'START_REVIEW.py'),
                               ('tree_review_view.js', 'tree_review_view.js'), ('tree_review.html', 'TREE_REVIEW.html')]:
            self.assertEqual((self.source / source).read_bytes(), (self.root / target).read_bytes())
        self.assertTrue((self.root / 'tree_profile.tsv').exists())
        self.assertTrue((self.root / 'tree_placement_taxonomy.tsv').exists())
        self.assertEqual(len(json.loads((self.root / 'iqtree_manifest.json').read_text())), len(IQTREE_DOWNLOADS))

    def test_bundle_launcher_imports_without_package_install(self):
        result = subprocess.run([sys.executable, 'START_REVIEW.py', '--help'], cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr); self.assertIn('--no-browser', result.stdout)

    def test_persistence_survives_restart_and_never_migrates_v1(self):
        legacy = self.root / '.review' / 'tree_review_v1.json'; legacy.parent.mkdir(); legacy.write_text('{"manual":15}')
        state = ReviewState(self.root); review = default_review(self.manifest)
        request = review['panels']['test']['groups']['family']
        preview = validate_panel(self.manifest, review['panels']['test'])
        request.update(approved=True, confirmed_signature=preview['groups']['family']['signature'])
        state.save(review)
        restarted = ReviewState(self.root)
        self.assertEqual(prepare_overlay(restarted.reviews(), self.manifest)['panels']['test']['models']['model_A']['points'], 12)
        self.assertEqual(legacy.read_text(), '{"manual":15}')
        self.assertTrue(list((self.root / '.review' / 'history').glob('*.json')))

    def test_seed_roundtrip_revalidates_before_restore(self):
        review = default_review(self.manifest); item = review['panels']['test']
        group = item['groups']['family']; group.update(approved=True,
            confirmed_signature=validate_panel(self.manifest, item)['groups']['family']['signature'])
        with tempfile.TemporaryDirectory() as output:
            manifest = write_bundle(output, panels(), initial_reviews=review)
            state = ReviewState(output)
            self.assertEqual(prepare_overlay(state.reviews(), manifest)['panels']['test']['models']['model_A']['points'], 12)
        changed = panels(); changed['test']['groups']['family']['species_tree'] = '((A,C),(B,D),E);'
        with tempfile.TemporaryDirectory() as output:
            manifest = write_bundle(output, changed, initial_reviews=review)
            self.assertIsNone(prepare_overlay(ReviewState(output).reviews(), manifest)['panels']['test']['models']['model_A']['points'])

    def test_rebuild_uses_selected_rows_and_keeps_provenance(self):
        item = {'panel': 'test', 'seal': 'sealed-source', 'group': 'family', 'inputSha256': 'stable-input',
            'selected_ids': list('abcd'), 'mask': {'mode': 'all'}}
        def runner(argv, **kwargs):
            self.assertIn('-keep-ident', argv); self.assertNotIn('-o', argv)
            (Path(kwargs['cwd']) / 'ml.treefile').write_text('((a,b)90,(c,d)90);')
            return SimpleNamespace(returncode=0, stdout='Success', stderr='')
        with patch('genecontext.tree_review_service._iqtree', return_value=Path('/fake/iqtree')), \
             patch('genecontext.tree_review_service.subprocess.run', return_value=SimpleNamespace(stdout='IQ-TREE test')):
            result = rebuild(self.root, self.manifest, item, runner=runner)
        self.assertFalse(result['approved']); self.assertEqual(result['source']['provenance']['selected_ids'], list('abcd'))
        req = default_review(self.manifest)['panels']['test']; req['groups']['family'].update(result)
        self.assertFalse(validate_panel(self.manifest, req)['groups']['family']['requires_rebuild'])

    def test_ui_does_not_implement_scoring_or_manual_split_controls(self):
        js = (self.root / 'tree_review_view.js').read_text(); html = (self.root / 'TREE_REVIEW.html').read_text()
        self.assertIn('/api/assess', js); self.assertIn('confirmed_signature', js)
        self.assertNotIn('Candidate-side', html); self.assertNotIn('modelSplits', js)
        self.assertIn('explicitly unrooted', html); self.assertIn('Rebuild tree', html)


class LocalHTTP(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name); write_bundle(self.root, panels())
        self.server = create_server(self.root); self.thread = threading.Thread(target=self.server.serve_forever, daemon=True); self.thread.start()

    def tearDown(self):
        self.server.shutdown(); self.server.server_close(); self.thread.join(); self.tmp.cleanup()

    def call(self, path, body=None, headers=None):
        connection = HTTPConnection('127.0.0.1', self.server.server_port)
        data = None if body is None else json.dumps(body)
        opts = dict(headers or {})
        if data is not None: opts.setdefault('Content-Type', 'application/json')
        connection.request('GET' if body is None else 'POST', path, data, opts)
        response = connection.getresponse(); content = response.read(); status = response.status; connection.close()
        try: content = json.loads(content)
        except ValueError: pass
        return status, content

    def test_same_origin_and_token_required_for_mutation(self):
        code, init = self.call('/api/manifest'); self.assertEqual(code, 200)
        item = default_review(init['manifest'])['panels']['test']
        self.assertEqual(self.call('/api/save', item)[0], 400)
        self.assertEqual(self.call('/api/save', item, {'X-GCF-Token': init['token'], 'Origin': 'https://evil.example'})[0], 400)
        self.assertEqual(self.call('/api/save', item, {'X-GCF-Token': init['token']})[0], 200)

    def test_host_guard_and_hidden_files(self):
        self.assertEqual(self.call('/api/manifest', headers={'Host': 'evil.example'})[0], 403)
        for path in ['/.review/tree_review_v2.json', '/%2e%2e/secret', '/START_REVIEW.py', '/tree_inputs.json']:
            self.assertNotEqual(self.call(path)[0], 200)
        self.assertEqual(self.call('/TREE_REVIEW.html')[0], 200)
        self.assertEqual(self.call('/tree_review_view.js')[0], 200)

    def test_linked_stage4_assets_and_same_origin_frame_are_available(self):
        (self.root / 'report.html').write_text('<html><script>window.ready=true</script><iframe src="TREE_REVIEW.html"></iframe></html>')
        (self.root / 'trees').mkdir(); (self.root / 'trees' / 'tree.nwk').write_text('(a,b);')
        (self.root / 'vendor').mkdir(); (self.root / 'vendor' / 'viewer.js').write_text('window.viewer=true')
        for path in ['/report.html', '/trees/tree.nwk', '/vendor/viewer.js']:
            self.assertEqual(self.call(path)[0], 200)
        with tempfile.TemporaryDirectory() as outside:
            secret = Path(outside) / 'secret.txt'; secret.write_text('private')
            (self.root / 'linked.txt').symlink_to(secret)
            self.assertNotEqual(self.call('/linked.txt')[0], 200)

    def test_post_preview_export_roundtrip_and_forged_values(self):
        _, init = self.call('/api/manifest'); headers = {'X-GCF-Token': init['token']}
        item = default_review(init['manifest'])['panels']['test']
        _, preview = self.call('/api/preview', item, headers)
        item['groups']['family'].update(approved=True, confirmed_signature=preview['groups']['family']['signature'], points=15)
        code, saved = self.call('/api/export', item, headers)
        self.assertEqual(code, 200); self.assertEqual(saved['overlay']['panels']['test']['models']['model_A']['points'], 12)
        _, restored = self.call('/api/state')
        self.assertEqual(restored['overlay']['panels']['test']['models']['model_A']['points'], 12)
        self.assertNotIn('points', restored['reviews']['panels']['test']['groups']['family'])

    def test_malformed_json_object_does_not_crash_handler(self):
        _, init = self.call('/api/manifest')
        self.assertEqual(self.call('/api/assess', [], {'X-GCF-Token': init['token']})[0], 400)

    def test_real_browser_review_persistence_and_invalidation(self):
        (self.root / 'report.html').write_text('<!doctype html><title>Stage4</title><script>window.reportReady=true</script><iframe src="TREE_REVIEW.html?panel=test&group=family"></iframe>')
        node = shutil.which('node') or str(Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node')
        if not Path(node).is_file(): self.skipTest('Node is unavailable')
        project = Path(__file__).resolve().parents[1]
        probe = subprocess.run([node, '-e', "require('playwright')"], cwd=project, capture_output=True, text=True)
        if probe.returncode: self.skipTest('Existing Playwright installation unavailable')
        script = r'''
const assert = require('assert/strict');
const {chromium} = require('playwright');
(async()=>{
 const browser = await chromium.launch({headless:true});
 try {
  const page=await browser.newPage(); const errors=[]; page.on('pageerror',e=>errors.push(e.message));
  await page.goto(process.argv[1]+'?panel=test&group=family');
  const ready=()=>page.waitForFunction(()=>!document.getElementById('approved').disabled);
  const included=()=>page.locator('#models tr').first().locator('td').nth(2);
  await ready(); assert.equal(await page.locator('#models tr').count(),5);
  assert.equal(await included().textContent(),'NA');
  await page.locator('#approved').check(); await page.waitForFunction(()=>document.querySelector('#models tr td:nth-child(3)').textContent==='12.000');
  await page.locator('#save').click(); await page.waitForFunction(()=>document.getElementById('activity').textContent==='Saved');
  await page.reload(); await ready(); assert.equal(await page.locator('#approved').isChecked(),true);
  await page.locator('#reviewOutgroup').selectOption(['a']); await page.locator('#reviewBasis').fill('Declared outgroup');
  await ready(); assert.equal(await page.locator('#approved').isChecked(),false); assert.equal(await included().textContent(),'NA');
  assert.equal(await page.locator('#models tr').first().locator('td').nth(1).textContent(),'12.000');
  await page.locator('#approved').check(); await page.waitForFunction(()=>document.querySelector('#models tr td:nth-child(3)').textContent==='12.000');
  await page.locator('#save').click(); await page.waitForFunction(()=>document.getElementById('activity').textContent==='Saved');
  await page.reload(); await ready(); assert.equal(await page.locator('#approved').isChecked(),true);
  assert.deepEqual(await page.locator('#reviewOutgroup').evaluate(e=>Array.from(e.selectedOptions,o=>o.value)),['a']);
  await page.locator('#stage3 summary').click();
  await page.locator('#members input[value="b"]').uncheck();
  await page.waitForFunction(()=>document.getElementById('gate').textContent.includes('rebuild first'));
  assert.equal(await page.locator('#approved').isDisabled(),true); assert.equal(await included().textContent(),'NA');
  await page.locator('#source').selectOption('other'); await ready();
  assert.equal(await page.locator('#models tr').first().locator('td').nth(1).textContent(),'0.000');
  assert.equal(await included().textContent(),'NA');
  await page.locator('#approved').check(); await page.waitForFunction(()=>document.querySelector('#models tr td:nth-child(3)').textContent==='0.000');
  await page.locator('#import').setInputFiles({name:'legacy.json',mimeType:'application/json',buffer:Buffer.from('{"schema":"gcf_tree_review_v1"}')});
  await page.waitForFunction(()=>document.getElementById('notice').textContent.includes('legacy manual T'));
  await page.goto(new URL('report.html',process.argv[1]).href); await page.waitForFunction(()=>window.reportReady===true);
  await page.frameLocator('iframe').locator('h1').waitFor();
  assert.equal(await page.frameLocator('iframe').locator('h1').textContent(),'Tree review');
  assert.deepEqual(errors,[]);
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
'''
        result = subprocess.run([node, '-e', script, 'http://127.0.0.1:' + str(self.server.server_port) + '/TREE_REVIEW.html'],
            cwd=project, capture_output=True, text=True, timeout=60)
        if 'Executable doesn\'t exist' in result.stderr: self.skipTest('Existing Chromium binary unavailable')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__': unittest.main()
