/* Several physical copies, one genome and one gene-comparison group. */
const fs=require('node:fs'),assert=require('node:assert/strict'),{JSDOM,VirtualConsole}=require('jsdom');
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(process.argv[2],'utf8'),{runScripts:'dangerously',url:'http://isolated-copies.invalid/',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}});
const w=dom.window,d=w.document,run=s=>w.eval(s),before=run('JSON.stringify(DATA)');
const change=(node,value)=>{node.value=value;node.dispatchEvent(new w.Event('change'))},control=name=>d.querySelector(`[aria-label="${name}"]`);
for(const r of run('DATA.regions')){change(control('Gene decision '+r.region_id),'unresolved');control('Save gene '+r.region_id).click()}
const rid=run(`DATA.regions.find(r=>r.genome_id.startsWith('Alligator')&&mainCandidates(r)[0]?.protein_length===310).region_id`);
const primary=control('Main model '+rid).value;
change(control('Gene decision '+rid),'inspect');control('Save gene '+rid).click();
const candidate=run(`rankedModels(own(regions.get(${JSON.stringify(rid)}))).find(m=>!m.overlaps_anchor&&m.coding_status==='complete_predicted'&&!m.window_edge&&m.target_support.length&&m.identification.proposed_name==='NANOG').model_id`);
control('Add gene to '+rid).click();
const addedOption=[...control('Gene to add').options].find(o=>run(`modelConstraint(models.get(${JSON.stringify(o.value)})).physical_id`)===run(`modelConstraint(models.get(${JSON.stringify(candidate)})).physical_id`));
assert.ok(addedOption);change(control('Gene to add'),addedOption.value);change(control('Added gene purpose'),'copy');
[...d.querySelectorAll('#addGenePanel button')].find(b=>b.textContent==='Add to table').click();
let extra=addedOption.value;assert.equal(control('Gene purpose '+extra).value,'copy');
assert.equal(control('Review comparison group '+extra).value,'NANOG');
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect'&&c.role==='target').length"),1);
change(control('Gene decision '+extra),'inspect');control('Save gene '+extra).click();
assert.equal(run('validate().length'),0,'A supported off-seed paralog needs no fake family split or expert exception');
assert.equal(run(`state.models[${JSON.stringify(extra)}].selection_kind`),'additional_copy');
assert.equal(run(`state.models[${JSON.stringify(extra)}].target_override`),false);
assert.equal(run("Object.values(state.models).filter(c=>c.decision==='inspect'&&c.role==='target'&&c.comparison_group==='NANOG').length"),2);
assert.notEqual(run(`copyLabel(models.get(${JSON.stringify(primary)}))`),run(`copyLabel(models.get(${JSON.stringify(extra)}))`));
const good=run('JSON.parse(JSON.stringify(state))');
for(const key of [rid,extra]){const n=control('Copy label '+key);n.value='same-copy';n.dispatchEvent(new w.Event('input'));control('Save gene '+key).click()}
assert.ok(run('validate().some(x=>x.includes("distinct copy labels"))'));
w.__good=good;run('load(__good);renderReviewQueue()');
if(process.argv[3])fs.writeFileSync(process.argv[3],run('JSON.stringify(state,null,2)')+'\n');
// Selecting an observation of the same physical copy in another reviewed
// window is still a duplicate, even with an unrelated display label.
const duplicate=run(`DATA.models.find(m=>m.region_id!==${JSON.stringify(rid)}&&modelConstraint(m).physical_id===modelConstraint(models.get(${JSON.stringify(extra)})).physical_id&&m.overlaps_anchor)`);
assert.ok(duplicate);change(control('Gene decision '+duplicate.region_id),'inspect');control('Save gene '+duplicate.region_id).click();
assert.ok(run('validate().some(x=>x.includes("more than once"))'));
assert.equal(run('JSON.stringify(DATA)'),before);assert.deepEqual(errors,[]);
console.log(JSON.stringify({status:'passed',two_copies_one_genome_one_group:true,added_copy_not_preapproved:true,off_seed_supported_copy_without_override:true,copy_labels_distinct:true,duplicate_labels_blocked:true,repeated_physical_copy_still_blocked:true,evidence_unchanged:true,dom_errors:errors},null,2));
dom.window.close();
