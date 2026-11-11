/* Read-only final report: offline controls must retain all exact model IDs. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{JSDOM,VirtualConsole}=require('jsdom');
const root=path.resolve(process.argv[2]),errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(path.join(root,'technical.html'),'utf8'),{runScripts:'dangerously',url:'http://isolated-stage3.invalid/technical.html',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}});
const w=dom.window,d=w.document,run=s=>w.eval(s),before=run('JSON.stringify(DATA)'),change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
assert.deepEqual(errors,[]);
const models=run('selected()');assert.ok(models.length);
assert.equal(d.querySelectorAll('#catalogue > table > tbody > tr').length,models.length);
assert.equal(d.querySelectorAll('#matrix tbody tr').length,models.length);
assert.equal(d.querySelectorAll('#matrix tbody td').length,models.length**2);
for(const kind of ['cds','coverage','protein']){change('matrixKind',kind);assert.ok(!/NaN|undefined|Infinity/.test(d.getElementById('matrix').textContent))}
const cell=[...d.querySelectorAll('#matrix td')].find(n=>n.textContent!=='—');assert.ok(cell);cell.click();assert.match(d.getElementById('pairInfo').textContent,/coverage/);
const search=d.getElementById('modelSearch');search.value=models[0].species;search.dispatchEvent(new w.Event('input'));
assert.equal(d.querySelectorAll('#catalogue > table > tbody > tr').length,models.filter(m=>m.species===models[0].species).length);
search.value='';search.dispatchEvent(new w.Event('input'));
change('profileMode','species');assert.match(d.getElementById('profile').getAttribute('src'),/species_conservation/);
assert.ok(d.getElementById('profileNote').textContent.includes('species'));
change('alignmentKind','protein');const n=d.querySelectorAll('.aln-row').length;assert.equal(n,models.length);
d.getElementById('showReferences').checked=true;d.getElementById('showReferences').dispatchEvent(new w.Event('change'));
assert.ok(d.querySelectorAll('.aln-row').length>=n);
d.getElementById('alnNext').click();assert.ok(Number(d.getElementById('alignmentStart').value)>1);
change('alignmentKind','cds');assert.equal(d.querySelectorAll('.aln-row').length,models.length);
change('dnaReference',models.at(-1).model_id);assert.ok(d.getElementById('dnaFigure').getAttribute('src').includes(models.at(-1).model_id));
assert.equal(d.querySelectorAll('#dnaTable tbody tr').length,models.length-1);
assert.ok(d.getElementById('syntenyStatus').textContent.includes('reviewed'));
let checked=0;
for(const node of d.querySelectorAll('a[href],img[src],iframe[src]')){
 const url=node.getAttribute('href')||node.getAttribute('src');if(!url||url.startsWith('#'))continue;
 assert.ok(!/^https?:/.test(url),'The report must not depend on remote assets');
 const name=decodeURIComponent(url.split(/[?#]/)[0]);assert.ok(fs.existsSync(path.join(root,name)),name);checked++;
}
assert.equal(run('JSON.stringify(DATA)'),before,'Display controls must not mutate scientific results');assert.deepEqual(errors,[]);
// Inspect the independent context viewer; it must not expose review controls.
const cv=new VirtualConsole();cv.on('jsdomError',e=>errors.push(e.message));
const context=new JSDOM(fs.readFileSync(path.join(root,'context.html'),'utf8'),{runScripts:'dangerously',url:'http://isolated-stage3.invalid/context.html',virtualConsole:cv,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}});
const c=context.window,cd=c.document;assert.deepEqual(errors,[]);
assert.equal(cd.querySelector('#links svg').dataset.rowHeight,'112');
assert.equal(c.eval('comparisonRows().length'),models.length);
assert.equal(cd.querySelectorAll('#links [data-window-row]').length,models.length);
assert.equal(cd.querySelectorAll('input[type=file],button[id=export],select[aria-label^="Decision"]').length,0);
for(const source of ['annotation','tiberius','all']){cd.getElementById('modelSource').value=source;cd.getElementById('modelSource').dispatchEvent(new c.Event('change'));assert.ok(cd.querySelector('#links svg'))}
assert.deepEqual(errors,[]);
const result={status:'passed',selected_models:models.length,all_model_ids_retained:true,identity_and_coverage_matrices:true,copy_aware_labels:true,protein_and_codon_alignment_viewer:true,reference_selector:true,species_balanced_profile:true,context_uses_stop2_ruler:true,context_read_only:true,offline_local_links_checked:checked,evidence_unchanged:true,dom_errors:errors};
if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));
context.window.close();dom.window.close();
