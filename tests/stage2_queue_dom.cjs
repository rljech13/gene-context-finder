/* STOP-2 visual-review handoff. Decisions exist only in isolated test memory. */
const fs=require('node:fs'),assert=require('node:assert/strict'),{JSDOM,VirtualConsole}=require('jsdom');
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(process.argv[2],'utf8'),{runScripts:'dangerously',url:'http://isolated-review.invalid/',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}});
const w=dom.window,d=w.document,run=s=>w.eval(s),before=run('JSON.stringify(DATA)');
const change=(node,value)=>{node.value=value;node.dispatchEvent(new w.Event('change'))};
const control=name=>d.querySelector(`[aria-label="${name}"]`);
assert.deepEqual(errors,[]);
assert.equal(d.querySelectorAll('#reviewQueue tr[data-primary="true"]').length,run('DATA.regions.length'));
assert.equal(d.querySelectorAll('#reviewQueue tr[data-primary="false"]').length,0);
assert.equal(run('Object.keys(state.models).length'),0,'Model proposals must never approve themselves');
assert.equal(run('Object.keys(state.regions).length'),0);
assert.equal(d.getElementById('other').value,'all');assert.equal(d.getElementById('modelLibrary').open,false);
// Save every other locus unresolved using the new table, without navigating it.
for(const r of run('DATA.regions')){change(control('Gene decision '+r.region_id),'unresolved');control('Save gene '+r.region_id).click()}
assert.equal(run('validate().length'),0,'Visual unresolved decisions produce a valid empty handoff');
const rid=run("DATA.regions.find(r=>{const m=models.get(reviewPlan(r).primary_model);return m&&m.coding_status==='complete_predicted'&&!m.window_edge&&!modelConstraint(m).target_warnings.length&&own(r).some(n=>modelConstraint(n).physical_id!==modelConstraint(m).physical_id&&n.coding_status==='complete_predicted'&&!n.window_edge)}).region_id");
const primary=control('Main model '+rid).value;
change(control('Gene decision '+rid),'inspect');control('Save gene '+rid).click();
assert.equal(run('validate().length'),0);
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect').length"),1);
// Adding one neighboring gene must not add any of the other window models.
control('Add gene to '+rid).click();const extra=control('Gene to add').value;assert.ok(extra);
const usable=run(`rankedModels(own(regions.get(${JSON.stringify(rid)}))).find(m=>m.coding_status==='complete_predicted'&&!m.window_edge&&modelConstraint(m).physical_id!==modelConstraint(models.get(${JSON.stringify(primary)})).physical_id).model_id`);
const option=[...control('Gene to add').options].find(o=>run(`modelConstraint(models.get(${JSON.stringify(o.value)})).physical_id`)===run(`modelConstraint(models.get(${JSON.stringify(usable)})).physical_id`));
change(control('Gene to add'),option.value);const added=option.value;
[...d.querySelectorAll('#addGenePanel button')].find(b=>b.textContent==='Add to table').click();
assert.equal(d.querySelectorAll('#reviewQueue tr[data-primary="false"]').length,1);
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect').length"),1,'An added proposal remains unapproved');
change(control('Added model '+added),usable);
change(control('Gene decision '+usable),'inspect');control('Save gene '+usable).click();
assert.equal(run('validate().length'),0);
assert.equal(run(`state.models[${JSON.stringify(usable)}].role`),'neighbor');
assert.equal(control('Main model '+rid).value,primary,'Adding a neighbor preserves the primary CDS');
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect'&&c.role==='target').length"),1);
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect'&&c.role==='neighbor').length"),1);
const pendingNote=control('Review comment '+rid);pendingNote.value='Isolated control: unsaved visual-review note';pendingNote.dispatchEvent(new w.Event('input'));
assert.equal(run(`state.regions[${JSON.stringify(rid)}].saved`),false,'An edited saved row must block export until saved again');
run('load(JSON.parse(JSON.stringify(state)));renderReviewQueue()');
assert.equal(control('Review comment '+rid).value,'Isolated control: unsaved visual-review note','Pending form edits must survive navigation and draft reload');
control('Save gene '+rid).click();assert.equal(run('validate().length'),0);
// Source/model changes in a gene never approve a replacement implicitly.
const saved=run('JSON.parse(JSON.stringify(state))');
if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(saved,null,2)+'\n');
run('load(JSON.parse(JSON.stringify(state)));renderReviewQueue()');
assert.equal(d.querySelectorAll('#reviewQueue tr[data-primary="false"]').length,1);
assert.equal(control('Main model '+rid).value,primary);
assert.throws(()=>run(`load({...state,review_plan:{[${JSON.stringify(rid)}]:{primary_model:'fake',additional:[]}}})`),/another locus/);
const alternatives=[...control('Main model '+rid).options].filter(o=>o.value!==primary);
if(alternatives.length){
 change(control('Main model '+rid),alternatives[0].value);
 assert.equal(run(`state.models[${JSON.stringify(primary)}].decision`),'not_selected');
 assert.ok(!run(`state.models[${JSON.stringify(alternatives[0].value)}]?.decision==='inspect'`),'A replacement must be confirmed again');
 run('load(JSON.parse(JSON.stringify(state)))');
 w.__saved=saved;run('load(__saved);renderReviewQueue()');
}
run('window.__downloads=[];download=(name,text)=>window.__downloads.push({name,text})');d.getElementById('export').click();
assert.equal(run('window.__downloads.length'),1);assert.equal(run('JSON.stringify(DATA)'),before);
// Explicitly request sequence comparison of an added gene, with a separate
// label and documented exception. It must preserve the original target.
change(control('Gene purpose '+usable),'target');
control('Review comparison group '+usable).value='ADDED_GENE_CONTROL';
control('Review comment '+usable).value='Isolated software control: intentionally compare this additional off-anchor gene';
control('Review expert override '+usable).checked=true;
control('Save gene '+usable).click();
assert.equal(run('validate().length'),0);
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect'&&c.role==='target').length"),2);
if(process.argv[4])fs.writeFileSync(process.argv[4],run('JSON.stringify(state,null,2)')+'\n');
// An exclusion clears selected models and selected pairs, not merely the table.
change(control('Gene decision '+rid),'exclude');control('Save gene '+rid).click();
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect').length"),0);assert.equal(run('validate().length'),0);
assert.deepEqual(errors,[]);
run('load(JSON.parse(JSON.stringify(state)))');
console.log(JSON.stringify({status:'passed',main_rows:run('DATA.regions.length'),default_main_only:true,no_automatic_approvals:true,decisions_without_locus_switching:true,exact_main_and_added_neighbor_export:true,explicit_additional_gene_comparison:true,main_preserved_on_add:true,replacement_requires_confirmation:true,pending_edits_preserved_and_export_blocked:true,draft_roundtrip:true,foreign_model_plan_rejected:true,exclude_clears_selections:true,scientific_evidence_unchanged:true,dom_errors:errors},null,2));
dom.window.close();
