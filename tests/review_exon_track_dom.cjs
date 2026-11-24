const assert = require('node:assert/strict');
const fs = require('node:fs');
const {webcrypto} = require('node:crypto');
const {JSDOM} = require('jsdom');

const seed = {stage: 2, review_id: 'splice-test', visible: ['plus', 'minus', 'query'], sequences: [
  {id: 'plus', group: 'g', label: 'Plus strand', protein: 'MKP', cds: 'ATGAAACCCTAA',
    metadata: {strand: '+', cds_blocks: [[100, 102, 0], [102, 104, 1], [204, 209, 2]]}},
  {id: 'minus', group: 'g', label: 'Minus strand', protein: 'MKP', cds: 'ATGAAACCC',
    metadata: {strand: '-', cds_blocks: [[30, 35, 2], [200, 205, 1]]}},
  {id: 'query', group: 'g', label: 'Protein query', protein: 'MIKP', cds: '', reference: true},
]};
const alignment = {id: 'splice-alignment', phase: 'alignment', group: 'g', review_id: seed.review_id, ids: seed.visible,
  protein: {plus: 'M-KP', minus: 'M-KP', query: 'MIKP'},
  cds: {plus: 'ATG---AAACCC', minus: 'ATG---AAACCC'}, mafft_version: 'fixture'};
const html = '<meta name="gcf-review-flow" content="sequential-v2"><section id="interactiveReview"></section>' +
  '<script id="interactiveReviewData" type="application/json">' + JSON.stringify(seed) + '</script>';
const dom = new JSDOM(html, {url: 'http://127.0.0.1:8795/fixture.html', runScripts: 'outside-only'});
const w = dom.window, d = w.document;
Object.defineProperty(w, 'crypto', {value: webcrypto});
w.TextEncoder = TextEncoder;
w.fetch = async () => ({ok: true, json: async () => ({review_id: seed.review_id})});
w.localStorage.setItem('gcf-realign-selection-splice-test-saved-alignment-g', JSON.stringify(alignment));
w.eval(fs.readFileSync(require('node:path').join(__dirname, '../genecontext/review_realign.txt'), 'utf8'));
const rows = () => [...d.querySelectorAll('.interactive-aln-row code')].map(row => [...row.children]);
const legend = d.querySelector('.gcf-exon-legend');
assert.match(legend.textContent, /2 \/ 3 rows mapped/);
let [plus, minus, query] = rows();
assert.equal(plus.map(cell => cell.textContent).join(''), 'M-KP');
assert.equal(plus[1].classList.contains('gcf-exon-cell'), false);
assert.equal(plus[1].title, '');
assert.match(plus[2].title, /CDS exon 1 · CDS exon 2 · Intron 1 · 100 bp · after base 1/);
assert.match(plus[2].style.getPropertyValue('--splice-ticks'), /33\.333/);
assert.match(minus[2].title, /Intron 1 · 165 bp · after base 1/);
assert.equal(query.some(cell => cell.classList.contains('gcf-exon-cell')), false);
const kind = d.querySelector('[aria-label="Interactive alignment kind"]');
kind.value = 'cds';
kind.dispatchEvent(new w.Event('change'));
[plus, minus] = rows();
assert.equal(plus.map(cell => cell.textContent).join(''), 'ATGAAACCC');
assert.match(plus[4].title, /CDS exon 2 · Intron 1 · 100 bp$/);
assert.match(minus[4].title, /CDS exon 2 · Intron 1 · 165 bp$/);
assert.match(legend.textContent, /2 \/ 2 rows mapped/);
legend.querySelector('input').click();
assert.equal(d.querySelectorAll('.gcf-exon-cell').length, 0);
legend.querySelector('input').click();
assert.equal(d.querySelectorAll('.gcf-exon-cell').length, 18);
dom.window.close();
console.log('Exact plus/minus splice sites, split codon, gap columns, terminal stop and display toggle verified.');
