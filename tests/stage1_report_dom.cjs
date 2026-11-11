/* Optional DOM regression runner. No browser, network, or scientific tools.
 * npm install --prefix /tmp/gcf-ui jsdom
 * NODE_PATH=/tmp/gcf-ui/node_modules node tests/stage1_report_dom.cjs REPORT \
 *   [SOURCE_EVIDENCE] [LEGACY_DRAFT] [DIAGNOSTICS_JSON]
 */
'use strict';
const fs=require('node:fs');
const assert=require('node:assert/strict');
const {JSDOM,VirtualConsole}=require('jsdom');
const [file,sourceFile,draftFile,output]=process.argv.slice(2);
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const started=Date.now();
const dom=new JSDOM(fs.readFileSync(file,'utf8'),{runScripts:'dangerously',url:'http://localhost:8780/',virtualConsole:vc});
const w=dom.window,d=w.document,run=s=>w.eval(s),plain=x=>JSON.parse(JSON.stringify(x));
const loadMs=Date.now()-started;
assert.equal(errors.length,0,errors.join('\n'));
assert.ok(d.querySelectorAll('.candidate').length<=10);
assert.equal(run('draft().pairs.reduce((n,p)=>n+Object.keys(p.overrides).length,0)'),0);
assert.ok(!d.querySelector('#overview').textContent.includes('Most proposed'));
const data=run('DATA');
const before=JSON.stringify(data.loci);
if(sourceFile){
 const original=JSON.parse(fs.readFileSync(sourceFile,'utf8'));
 for(const key of ['loci','queries','genomes','filters','review_id','early_synteny'])assert.deepEqual(plain(data[key]),original[key],key+' was changed by transport');
}
let alternateProducts=0,unknownNames=0,displayCopies=0,syntenyViews=0;
for(const l of data.loci){
 if(!l.annotation_gene_name){assert.equal(run(`locusName(lmap.get(${JSON.stringify(l.locus_id)}))`),'Candidate '+l.locus_id);unknownNames++}
 for(const q of data.queries){
  w.__l=l;w.__qid=q.id;
  const anyPass=run('__l.chains.some(c=>c.query===__qid&&chainPasses(c))');
  const chosen=run('supportChain(__l,__qid)'),best=run('bestChain(__l,__qid)');
  if(anyPass){w.__c=chosen;assert.equal(run('chainPasses(__c)'),true);if(chosen.chain_id!==best.chain_id)alternateProducts++}
 }
}
for(const g of data.genomes){
 d.getElementById('genome').value=g.genome_id;run('shownLimit=10;render()');
 const cards=[...d.querySelectorAll('.candidate')];
 assert.ok(cards.length<=10,'Initial batch exceeded ten');
 // Closed HSP/CDS/ranking panels must not build their contents eagerly.
 for(const card of cards){assert.equal(card.querySelectorAll('.maps').length,0);assert.ok(!card.textContent.includes('Original bit-score discovery rank:'))}
 const own=run(`groups.get(${JSON.stringify(g.genome_id)})`);
 for(const l of own.slice(0,10)){
  if(!l.early_synteny?.comparisons?.length)continue;
  w.__l=l;const view=run('synteny(__l)');syntenyViews++;
  assert.equal(view.querySelectorAll('svg.neighbor-map').length,1,'Duplicate synteny diagram');
  assert.equal(view.querySelectorAll('select[aria-label="Overview genes shown"] option[value="20"]').length,0);
  for(const c of l.early_synteny.comparisons){const r=run(`refsById.get(${JSON.stringify(c.reference_id)})`);if(!r)continue;w.__r=r;
   const genes=run('displayGenes(__r,{...__r,neighbors_used:__l.early_synteny.neighbors_used},"10")');
   for(const side of ['upstream','downstream'])assert.ok(genes.filter(g=>!g.display_only&&g.target_side===side).length<=10);
  }
  const extras=view.querySelectorAll('polygon[data-match-type="display-only"]');displayCopies+=extras.length;
  // Every saved reference is represented exactly once, including controls.
  assert.equal(view.querySelectorAll('svg.neighbor-map').item(0).getAttribute('aria-label'),'Candidate neighborhood compared with all available query-source loci');
 }
}
assert.equal(JSON.stringify(data.loci),before,'Rendering mutated shared scientific objects');
let savedPairs=null,selectedLoci=null,sparseEdits=null,backendRows=null;
if(draftFile){
 const originalDraft=JSON.parse(fs.readFileSync(draftFile,'utf8'));w.__draft=originalDraft;run('loadDraft(__draft)');
 const migrated=plain(run('draft()'));
 savedPairs=migrated.pairs.filter(p=>p.committed).length;
 selectedLoci=migrated.pairs.reduce((n,p)=>n+p.selected.length,0);
 sparseEdits=migrated.pairs.reduce((n,p)=>n+Object.keys(p.overrides).length,0);
 for(const p of originalDraft.pairs){const loaded=migrated.pairs.find(x=>x.genome_id===p.genome_id);assert.deepEqual(loaded.selected,p.selected);assert.equal(loaded.action,p.action);assert.equal(loaded.reason,p.reason);assert.equal(loaded.notes,p.notes);assert.equal(loaded.committed,p.committed)}
 const rows=plain(run('exportRows()'));
 for(const r of rows){assert.ok(!r.notes.includes('loci reviewed'));if(r.selection_status==='not_selected_for_run')assert.ok(r.notes.includes('individual review is not recorded'));assert.equal(r.display_filters,'unknown_legacy')}
 // Persist a fresh global filter snapshot, round-trip it, and retain old unknown history.
 run("$('coverage').value=.35;$('onlyPassing').checked=true");const newDraft=plain(run('draft()'));
 run("$('coverage').value=.99;$('onlyPassing').checked=false");w.__draft=newDraft;run('loadDraft(__draft)');
 assert.equal(run('displayFilters().min_query_coverage'),.35);assert.equal(run('displayFilters().only_passing'),true);
 // A newly saved pair records its exact filters; TSV round-trips selected loci.
 const selectedPair=migrated.pairs.find(p=>p.action==='inspect'&&p.selected.length);
 if(selectedPair){d.getElementById('genome').value=selectedPair.genome_id;run('render()');d.getElementById('savePair').click();assert.equal(run('currentState().decision_filters.min_query_coverage'),.35)}
 // Only a complete saved review can be exported to TSV. Complete pending pairs
 // in this isolated test copy; never alter the user's portable draft.
 run("for(const s of state.values())if(!s.committed){s.action=s.selected.length?'inspect':'unresolved';s.reason=s.selected.length?'incomplete_evidence':'insufficient_evidence';s.committed=true;s.decision_filters=displayFilters()}");
 assert.equal(run('validateAll().length'),0);
 const beforeTSV=plain(run('draft()'));w.__rows=run('exportRows()');
 const tsv=run("[columns.join('\\t'),...__rows.map(r=>columns.map(k=>String(r[k]??'').replace(/[\\t\\r\\n]/g,' ')).join('\\t'))].join('\\n')+'\\n'");w.__tsv=tsv;run('loadTSV(__tsv)');
 const afterTSV=plain(run('draft()'));for(const p of beforeTSV.pairs){const loaded=afterTSV.pairs.find(x=>x.genome_id===p.genome_id);assert.deepEqual(loaded.selected,p.selected);assert.deepEqual(loaded.decision_filters,p.decision_filters||null)}
 backendRows=plain(run('exportRows()'));
 // The loaded user draft is only an in-memory test; no file or browser storage was changed.
}
assert.equal(errors.length,0,errors.join('\n'));
const result={dom_errors:errors,html_bytes:fs.statSync(file).size,jsdom_initialization_ms:loadMs,
 scientific_transport_equal:!!sourceFile,rendering_preserves_evidence:true,
 alternative_passing_products:alternateProducts,unnamed_candidates_correctly_labeled:unknownNames,
 tested_synteny_views:syntenyViews,display_only_gene_occurrences:displayCopies,
 migrated_saved_pairs:savedPairs,migrated_selected_loci:selectedLoci,migrated_sparse_boundary_edits:sparseEdits,
 legacy_draft_choices_preserved:!!draftFile,filter_snapshot_and_tsv_round_trip:!!draftFile};
if(output){fs.writeFileSync(output,JSON.stringify(result,null,2)+'\n');if(backendRows)fs.writeFileSync(output.replace(/\.json$/,'.synthetic_completed_review.rows.json'),JSON.stringify(backendRows,null,2)+'\n')}
console.log(JSON.stringify(result,null,2));dom.window.close();
