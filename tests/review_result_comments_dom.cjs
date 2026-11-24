const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {JSDOM} = require('jsdom');

const source = fs.readFileSync(path.join(__dirname, '../genecontext/alignment_review_steps.txt'), 'utf8');
const seed = {review_id: 'comment-test', source_review_id: 'native-review', sequences: []};
const prefix = 'gcf-realign-selection-' + seed.review_id;
const html = `<meta name="gcf-review-flow" content="sequential-v2">
  <section id="interactiveReview"><h2>Alignment</h2><p>Review</p>
    <select aria-label="Interactive gene group"><option>Test</option></select>
    <select aria-label="Interactive alignment kind"><option>protein</option></select>
    <div id="interactiveAlignment"></div><button id="realignSelected">Align</button>
    <p id="realignStatus"></p></section>
  <script id="interactiveReviewData" type="application/json">${JSON.stringify(seed)}</script>`;
const ids = ['first', 'second'];

function open(values = {}) {
  const dom = new JSDOM(html, {url: 'http://127.0.0.1:8795/report.html', runScripts: 'outside-only'});
  const w = dom.window;
  const d = w.document;
  const host = d.getElementById('interactiveReview');
  let alignment, result;
  for (const [key, value] of Object.entries(values)) w.localStorage.setItem(key, value);
  w.GCFAlignmentWorkflow = {
    current: () => ({ids, result}), savedAlignment: () => alignment,
    treeIDs: () => ids, alignmentMatchesModels: () => true,
  };
  w.eval(source);
  return {
    w, alignmentNote: d.getElementById('alignmentReviewComment'), treeNote: d.getElementById('treeReviewComment'),
    show(alignmentId, treeId) {
      alignment = alignmentId ? {id: alignmentId, group: 'Test', ids} : undefined;
      result = treeId ? {id: treeId, phase: 'tree', alignment_id: alignmentId, ids, phylogeny: {status: 'computed'}} : undefined;
      host.dispatchEvent(new w.Event('gcf-review-state'));
    },
    type(input, value) {
      input.value = value;
      input.dispatchEvent(new w.Event('input', {bubbles: true}));
    },
  };
}

const page = open();
assert.equal(page.alignmentNote.parentElement.hidden, true);
assert.equal(page.treeNote.parentElement.hidden, true);
page.show('alignment-1', 'tree-1');
assert.equal(page.alignmentNote.required, false);
assert.equal(page.treeNote.required, false);
page.type(page.alignmentNote, 'Inspect a gap\nCheck CDS');
page.type(page.treeNote, 'Review that clade');
page.show('alignment-1', 'tree-2');
assert.equal(page.alignmentNote.value, 'Inspect a gap\nCheck CDS');
assert.equal(page.treeNote.value, '');
page.show('alignment-2');
assert.equal(page.alignmentNote.value, '');
assert.equal(page.treeNote.parentElement.hidden, true);
assert.equal(page.w.localStorage.getItem(prefix + '-tree-comment-tree-1'), 'Review that clade');

const reload = open(Object.fromEntries(Object.keys(page.w.localStorage).map(key => [key, page.w.localStorage.getItem(key)])));
reload.show('alignment-1', 'tree-1');
assert.equal(reload.alignmentNote.value, 'Inspect a gap\nCheck CDS');
assert.equal(reload.treeNote.value, 'Review that clade');
assert.equal(reload.treeNote.parentElement.hidden, false);
reload.type(reload.treeNote, '');
assert.equal(reload.w.localStorage.getItem(prefix + '-tree-comment-tree-1'), '');
page.w.close();
reload.w.close();
console.log('Optional alignment and tree comments survive reload and follow their exact result IDs.');
