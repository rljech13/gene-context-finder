/* Optional local DOM test. NODE_PATH=/tmp/gcf_report_ui_test/node_modules node tests/stage2_report_dom.cjs REPORT [RESULT_JSON] */
const fs=require('node:fs'),assert=require('node:assert/strict');
const {JSDOM,VirtualConsole}=require('jsdom');
(async()=>{
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(process.argv[2],'utf8'),{runScripts:'dangerously',url:'http://localhost:8781/',virtualConsole:vc});
 const w=dom.window,d=w.document,run=s=>w.eval(s),tick=()=>new Promise(r=>setImmediate(r));await tick();
 assert.deepEqual(errors,[]);assert.equal(run('Object.keys(state.models).length'),0);assert.equal(run('Object.keys(state.regions).length'),0);
 const evidenceBefore=JSON.stringify(run('DATA'));let reviewed=0,selectedModels=0;
 // Display controls must neither move the focal gene away from the centre in
 // order view nor turn schematic spacing into an apparent base-pair scale.
 const dispatch=id=>d.getElementById(id).dispatchEvent(new w.Event('change'));
 assert.equal(d.getElementById('modelSource').value,'all');
 assert.equal(d.getElementById('scale').value,'distance','Real genomic ruler is the default');

 let overviewSpans=0;
 for(const r of run('DATA.regions')){
  d.getElementById('region').value=r.region_id;dispatch('region');await tick();
  const svg=d.querySelector('#links .cluster-overview'),px=Number(svg.dataset.pxPerBp),rows=run('comparisonRows()'),pitch=Number(svg.dataset.rowHeight),height=Number(svg.getAttribute('viewBox').split(' ')[3]);
  assert.equal(pitch,112);assert.equal(height,112*rows.length+68,'Overview height depends only on locus count, not isoform count');
  const centre=Number(svg.querySelector('[data-ruler-bp="0"]').getAttribute('x1')),shown=run('neighborhood(current())'),expected=[...shown.left,...shown.anchor,...shown.right].flatMap(g=>g.map(m=>m.model_id)).sort();
  const shapes=[...svg.querySelectorAll('.gene-arrow[data-row="0"] .cluster-model-span')];
  assert.deepEqual(shapes.flatMap(n=>n.dataset.modelIds.split(' ')).sort(),Array.from(expected),'Overview retains every original model ID');
  assert.equal(svg.querySelectorAll('.cluster-gene-label[data-row="0"]').length,shown.left.length+shown.anchor.length+shown.right.length);
  const sign=r.strand==='-'?-1:1,mid=(r.anchor_start+r.anchor_end)/2;
  for(const shape of shapes){
   const a=Number(shape.dataset.start),b=Number(shape.dataset.end),points=shape.getAttribute('points').split(/\s+/).map(p=>p.split(',').map(Number)),xs=points.map(p=>p[0]),x=centre+Math.min((a-mid)*sign,(b-mid)*sign)*px;
   assert.ok(Math.abs(Math.min(...xs)-x)<1e-8);assert.ok(Math.abs(Math.max(...xs)-Math.min(...xs)-(b-a)*px)<1e-8,'Overview spans use true bp lengths');overviewSpans++;
  }
  const labels=[...svg.querySelectorAll('.cluster-gene-label')];
  for(let i=0;i<labels.length;i++)for(let j=i+1;j<labels.length;j++){
   const a=labels[i].dataset,b=labels[j].dataset;assert.equal(Number(a.left)<Number(b.right)&&Number(b.left)<Number(a.right)&&Number(a.top)<Number(b.bottom)&&Number(b.top)<Number(a.bottom),false);
  }
 }
 for(const source of ['annotation','tiberius']){
  d.getElementById('modelSource').value=source;dispatch('modelSource');
  for(const shape of d.querySelectorAll('.cluster-model-span'))for(const id of shape.dataset.modelIds.split(' ')){w.__checkId=id;assert.equal(run('comparisonRows().flatMap(diagramModels).find(m=>m.model_id===__checkId).source'),source)}
  for(const label of d.querySelectorAll('.cluster-gene-label'))assert.equal(label.dataset.sourceBadge,source==='annotation'?'A':'T');
 }
 d.getElementById('modelSource').value='all';dispatch('modelSource');
 d.getElementById('scale').value='models';dispatch('scale');
 let rulerModels=0,rulerBlocks=0,rulerStructures=0,sharedStructures=0;
 for(const r of run('DATA.regions')){
  d.getElementById('region').value=r.region_id;dispatch('region');await tick();
  const svg=d.querySelector('#links .coordinate-ruler'),px=Number(svg.dataset.pxPerBp),rows=run('comparisonRows()');
  assert.ok(px>0);
  const spans=[...svg.querySelectorAll('[data-window-row] rect')];
  spans.forEach((rect,i)=>assert.ok(Math.abs(Number(rect.getAttribute('width'))/(rows[i].window_end-rows[i].window_start)-px)<1e-12));
  const labels=[...svg.querySelectorAll('.ruler-gene-label')];
  for(let i=0;i<labels.length;i++)for(let j=i+1;j<labels.length;j++){
   const a=labels[i].dataset,b=labels[j].dataset;
   const overlap=Number(a.left)<Number(b.right)&&Number(b.left)<Number(a.right)&&Number(a.top)<Number(b.bottom)&&Number(b.top)<Number(a.bottom);
   assert.equal(overlap,false,'Ruler labels must occupy non-overlapping boxes');
  }
  const drawn=[...svg.querySelectorAll('.transcript-model[data-row="0"]')],expected=run('neighborhood(current())');
  const ids=[...expected.left,...expected.anchor,...expected.right].flatMap(g=>g.map(m=>m.model_id)).sort();
  assert.deepEqual(drawn.flatMap(n=>n.dataset.modelIds.split(' ')).sort(),Array.from(ids),'Every supplied alternative remains represented exactly once');
  assert.equal(svg.querySelectorAll('.ruler-gene-label[data-row="0"]').length,expected.left.length+expected.anchor.length+expected.right.length,'One name per displayed gene');
  const centre=Number(svg.getAttribute('viewBox').split(' ')[2])/2;
  for(const model of drawn){
   const id=model.dataset.modelId;w.__checkId=id;
   const m=run('diagramModels(current()).find(m=>m.model_id===__checkId)'),sign=r.strand==='-'?-1:1,mid=(r.anchor_start+r.anchor_end)/2;
   const memberIds=model.dataset.modelIds.split(' '),badge=model.querySelector('.model-source-badge');
   const members=memberIds.map(id=>{w.__checkId=id;return run('diagramModels(current()).find(m=>m.model_id===__checkId)')});
   for(const v of members){
    for(const key of ['contig','strand','start','end','cds_sha256','context_only','biotype'])assert.equal(v[key],m[key],'Collapsed models must share complete native geometry and fingerprint');
    assert.deepEqual(JSON.parse(JSON.stringify(v.cds_blocks)),JSON.parse(JSON.stringify(m.cds_blocks)));
    assert.deepEqual(JSON.parse(JSON.stringify(v.stop_blocks||[])),JSON.parse(JSON.stringify(m.stop_blocks||[])));
    assert.equal(v.gene_group,m.gene_group);
   }
   const sources=new Set(members.map(v=>v.source)),expectedBadge=sources.size===2?'A+T':sources.has('annotation')?'A':'T';
   assert.equal(badge.childNodes[0].textContent,expectedBadge);
   if(sources.size===2)sharedStructures++;
   rulerStructures++;rulerModels+=members.length;
   const expectedBlocks=m.cds_blocks.filter(([a,b])=>Math.min(b,r.window_end)>Math.max(a,r.window_start));
   const blocks=[...model.querySelectorAll('.cds-block')];assert.equal(blocks.length,expectedBlocks.length);
   for(const block of blocks){
    const a=Number(block.dataset.start),b=Number(block.dataset.end),x=centre+Math.min((a-mid)*sign,(b-mid)*sign)*px;
    assert.ok(Math.abs(Number(block.getAttribute('x'))-x)<1e-8,'CDS x must use native base-pair coordinates');
    assert.ok(Math.abs(Number(block.getAttribute('width'))-(b-a)*px)<1e-8,'Tiny CDS blocks must not be stretched to a minimum width');rulerBlocks++;
   }
  }
 }
 const basePx=Number(d.querySelector('#links .coordinate-ruler').dataset.pxPerBp);
 d.getElementById('genomeZoom').value='2';dispatch('genomeZoom');
 assert.ok(Math.abs(Number(d.querySelector('#links .coordinate-ruler').dataset.pxPerBp)/basePx-2)<1e-12);
 d.getElementById('genomeZoom').value='1';dispatch('genomeZoom');
 d.getElementById('scale').value='order';dispatch('scale');

 const checkSource=source=>{
  for(const arrow of d.querySelectorAll('#links .gene-arrow'))for(const id of arrow.dataset.modelIds.split(' ')){
   w.__checkId=id;assert.equal(run('comparisonRows().flatMap(diagramModels).find(m=>m.model_id===__checkId).source'),source);
  }
  for(const p of run('diagramPairs()'))for(const k of ['model_a','model_b']){
   w.__checkId=p[k];assert.equal(run('models.get(__checkId).source'),source);
  }
 };
 const usedPhysical=new Set();
 for(const r of run('DATA.regions')){
  d.getElementById('region').value=r.region_id;dispatch('region');await tick();
  assert.equal(run('diagramModels(current()).length'),run('own(current()).length+contextFeatureModels(current()).length'));
  const shown=run('neighborhood(current())');
  assert.equal(d.querySelectorAll('#links .gene-arrow[data-row="0"]').length,shown.left.length+shown.anchor.length+shown.right.length);
  const caps=[...d.querySelectorAll('#links [data-window-edge][data-row="0"]')];
  assert.equal(caps.length,2,'Both physical display edges must remain visible');
  assert.deepEqual(caps.map(n=>Number(n.dataset.coordinate)),r.strand==='-'?[r.window_end,r.window_start]:[r.window_start,r.window_end]);
  const count=d.querySelector('#links .gene-count[data-row="0"]');
  assert.equal(Number(count.dataset.shown),d.querySelectorAll('#links .gene-arrow[data-row="0"]').length);
  const own=run('own(current()).concat(contextFeatureModels(current()))');
  for(const source of ['annotation','tiberius'])assert.equal(Number(count.dataset[source+'Genes']),new Set(own.filter(m=>m.source===source).map(m=>m.gene_group)).size);
  const targets=[...d.querySelectorAll('#links .gene-arrow[data-row="0"][data-match-type="target"]')];
  if(targets.length===1){const xs=[...targets[0].querySelectorAll('polygon')].flatMap(n=>n.getAttribute('points').split(/\s+/).map(pair=>Number(pair.split(',')[0])));assert.equal((Math.min(...xs)+Math.max(...xs))/2,630)}
  for(const node of d.querySelectorAll('#links .gene-arrow[data-row="0"]')){
   w.__group=node.dataset.geneGroup;
   const expected=run('diagramModels(current()).filter(m=>m.gene_group===__group)');
   const shapes=[...node.querySelectorAll('.source-arrow')];
   assert.equal(shapes.length,new Set(expected.map(m=>m.source)).size);
   assert.deepEqual(shapes.flatMap(n=>n.dataset.modelIds.split(' ')).sort(),Array.from(expected,m=>m.model_id).sort());
  }
 }
 d.getElementById('region').value=run('DATA.regions[0].region_id');dispatch('region');
 d.getElementById('modelSource').value='annotation';dispatch('modelSource');checkSource('annotation');
 d.getElementById('modelSource').value='all';dispatch('modelSource');assert.equal(run('diagramModels(current()).length'),run('own(current()).length+contextFeatureModels(current()).length'));
 d.getElementById('modelSource').value='tiberius';dispatch('modelSource');checkSource('tiberius');
 if([...d.getElementById('other').options].some(o=>o.value==='all')){d.getElementById('other').value='all';dispatch('other')}
 const arrows=[...d.querySelectorAll('#links .gene-arrow')];assert.ok(arrows.length);
 const rowCount=run('comparisonRows().length');assert.equal(new Set(arrows.map(n=>n.dataset.row)).size,rowCount);
 for(const p of d.querySelectorAll('#links .neighbor-link')){
  assert.ok(!/NaN|Infinity/.test(p.getAttribute('d')));assert.equal(Number(p.dataset.toRow)-Number(p.dataset.fromRow),1);
 }
 const linked=d.querySelector('#links .neighbor-link');
 if(linked){linked.dispatchEvent(new w.Event('mouseenter'));assert.ok(d.querySelector('#links svg.link-focus'));assert.ok(linked.classList.contains('link-active'));linked.dispatchEvent(new w.Event('mouseleave'));assert.ok(!d.querySelector('#links svg.link-focus'))}
 d.getElementById('scale').value='distance';dispatch('scale');assert.match(d.getElementById('links').textContent,/kb/);
 const spans=[...d.querySelectorAll('#links [data-window-row] rect')],rows=run('comparisonRows()');
 const ratios=spans.map((span,i)=>Number(span.getAttribute('width'))/(rows[i].window_end-rows[i].window_start));
 assert.ok(ratios.every(v=>Math.abs(v-ratios[0])<1e-12),'Distance-view fragment spans must share one bp scale');
 d.getElementById('scale').value='order';dispatch('scale');assert.match(d.getElementById('links').textContent,/equal spacing/);
 assert.ok(!d.getElementById('tracks').closest('details').open,'Detailed CDS plot should start collapsed');
 run('window.__downloads=[];download=(name,text)=>window.__downloads.push({name,text})');
 for(const r of run('DATA.regions')){
  d.getElementById('region').value=r.region_id;d.getElementById('region').dispatchEvent(new w.Event('change'));await tick();
  assert.ok(d.getElementById('regionInfo').textContent.includes(r.contig));
  assert.equal(d.querySelectorAll('#tracks svg').length,1);
  const candidates=run('own(current())');
  const eligible=candidates.filter(m=>m.overlaps_anchor&&m.coding_status!=='disrupted'&&!usedPhysical.has(run(`modelConstraint(models.get(${JSON.stringify(m.model_id)})).physical_id`)));
  const m=eligible.find(m=>m.coding_status==='complete_predicted')||eligible[0];
  if(m){w.__mid=m.model_id;const card=[...d.querySelectorAll('#models details.gene')].find(n=>n.dataset.geneGroup===m.gene_group);card.open=true;await tick();
   const choice=d.querySelector('select[aria-label="Decision '+m.model_id+'"]');assert.ok(choice);choice.value='target';choice.dispatchEvent(new w.Event('change'));
   const note=d.querySelector('input[aria-label="Comment '+m.model_id+'"]');note.value='Software control only; not a user biological decision';note.dispatchEvent(new w.Event('input'));selectedModels++;
   const override=d.querySelector('input[aria-label="Expert target override '+m.model_id+'"]');override.checked=true;override.dispatchEvent(new w.Event('change'));
   const group=d.querySelector('input[aria-label="Comparison group '+m.model_id+'"]');group.value='isolated_control_'+r.region_id;group.dispatchEvent(new w.Event('input'));
   usedPhysical.add(run('modelConstraint(models.get(__mid)).physical_id'));
  }
  d.getElementById('regionDecision').value=m?'inspect':'unresolved';d.getElementById('regionDecision').dispatchEvent(new w.Event('change'));
  d.getElementById('regionNotes').value='Isolated software control';d.getElementById('regionNotes').dispatchEvent(new w.Event('input'));
  d.getElementById('saveRegion').click();assert.equal(run(`state.regions[${JSON.stringify(r.region_id)}].saved`),true);reviewed++;
 }
 assert.equal(run('validate().length'),0);d.getElementById('export').click();
 let downloads=run('window.__downloads');assert.equal(downloads.length,1,'Primary export is one authoritative JSON');const exported=JSON.parse(downloads.find(x=>x.name==='stage2.decisions.json').text);
 d.getElementById('exportTables').click();downloads=run('window.__downloads');assert.equal(downloads.length,4);
 assert.equal(exported.review_id,run('DATA.review_id'));assert.ok(downloads.find(x=>x.name==='models.decisions.tsv').text.includes('not_selected')||run('DATA.models.length')===selectedModels);
 run('load(JSON.parse(JSON.stringify(state)))');assert.equal(run('validate().length'),0);
 assert.throws(()=>run("load({...state,review_id:'wrong'})"),/another/);
 assert.throws(()=>run("load({...state,regions:{fake:{decision:'inspect',saved:true,notes:'x'}}})"),/Unknown/);
 assert.throws(()=>run("load({...state,pairs:{fake:{decision:'inspect'}}})"),/Unknown proposed pair/);
 for(const m of run('DATA.models'))if(m.annotation_name&&m.identification.proposed_name&&m.annotation_name.toUpperCase()!==m.identification.proposed_name.toUpperCase()){w.__mid=m.model_id;const label=run('groupLabel(DATA.models.filter(x=>x.region_id===models.get(__mid).region_id&&x.gene_group===models.get(__mid).gene_group))');assert.ok(label.includes(m.annotation_name));assert.ok(label.includes(m.identification.proposed_name))}
 assert.equal(JSON.stringify(run('DATA')),evidenceBefore,'UI mutated scientific evidence');assert.deepEqual(errors,[]);
 const result={default_fixed_height_cluster_overview:true,overview_spans_checked:overviewSpans,cds_detail_view_retained:true,compact_single_name_per_gene:true,all_original_models_represented:rulerModels,distinct_structure_drawings:rulerStructures,shared_annotation_tiberius_structures:sharedStructures,ruler_cds_blocks_checked:rulerBlocks,nonoverlapping_full_labels:true,exact_cds_bp_geometry:true,zoom_preserves_scale:true,dom_errors:errors,regions_reviewed_in_isolated_control:reviewed,models_selected_in_control:selectedModels,no_automatic_selections:true,evidence_unchanged:true,schematic_combined_separate_source_arrows:true,source_specific_links:true,annotation_and_combined_views:true,decision_export_and_reload:true,wrong_checkpoint_rejected:true,target_centring_and_relative_scale:true,multirow_neighbor_links_and_hover:true,window_edge_coordinates_in_target_orientation:true,source_gene_counts_and_shown_counts:true,physical_window_lengths_on_shared_bp_scale:true,note:'These in-memory test decisions are not user-approved and were not written into the report or evidence.'};
 if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));dom.window.close();
})().catch(e=>{console.error(e);process.exitCode=1});
