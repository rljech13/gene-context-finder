/* Independent stacked viewers share evidence and one set of model decisions. */
const fs=require('node:fs'),assert=require('node:assert/strict');
function checkComparisons(report,{JSDOM,VirtualConsole}=require('jsdom')){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const html=fs.readFileSync(report,'utf8');
 const open=storage=>new JSDOM(html,{runScripts:'dangerously',url:'http://localhost/stacked/',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){};for(const [k,v]of Object.entries(storage||{}))w.localStorage.setItem(k,v)}});
 const dom=open(),w=dom.window,d=w.document,run=s=>w.eval(s);
 const before=run('JSON.stringify(DATA)'),decisions=run('JSON.stringify(exportState())');
 const pane=n=>d.querySelectorAll('.comparison-pane')[n],control=(n,key)=>pane(n).querySelector(`[data-control="${key}"]`);
 const change=(n,key,value)=>{const node=control(n,key);node.value=value;node.dispatchEvent(new w.Event('change'))};
 d.getElementById('addComparison').click();d.getElementById('addComparison').click();
 assert.equal(d.querySelectorAll('.comparison-pane').length,3);
 assert.equal(d.querySelectorAll('#reviewQueue').length,1,'One shared model decision table');
 if(run('DATA.regions[0].gene')==='KLF4'){
  const options=[...control(0,'region').options].map(n=>n.textContent);
  assert.match(options[0],/^KLF4 · Gallus gallus/);
  assert.match(options[1],/^KLF2 · Gallus gallus/,'Locus names must describe the actual gene, not repeat the query-panel name');
  assert.match(pane(0).querySelector('h2').textContent,/KLF4/);
  assert.ok([...control(0,'links').querySelectorAll('[data-locus-name]')].some(n=>n.textContent.startsWith('KLF2')),'Each cluster row names its focal gene');
 }
 const ids=[...d.querySelectorAll('[id]')].map(n=>n.id);assert.equal(new Set(ids).size,ids.length,'No duplicate DOM IDs');
 assert.equal(d.getElementById('comparisonPanels').children.length,3,'Panes stack as separate full-width sections');
 const first=control(0,'links').innerHTML,third=control(2,'links').innerHTML;
 const rid=run('DATA.regions[1].region_id');
 change(1,'region',rid);change(1,'other','custom');control(1,'clearRows').click();
 const wanted=control(1,'rowChoices').querySelector('input');wanted.checked=true;wanted.dispatchEvent(new w.Event('change'));
 assert.equal(control(1,'links').querySelectorAll('.cluster-row-ruler').length,2,'Custom subset contains the focal locus and one checked locus');
 change(1,'genomeZoom','2');change(1,'modelSource','tiberius');change(1,'neighborsShown','5');
 assert.equal(control(0,'links').innerHTML,first,'Editing the second pane must not redraw the first');
 assert.equal(control(2,'links').innerHTML,third,'Editing the second pane must not redraw the third');
 assert.equal(control(0,'genomeZoom').value,'1');assert.equal(control(2,'modelSource').value,'all');
 change(2,'other',rid);change(2,'scale','models');
 assert.ok(control(2,'links').querySelector('.coordinate-ruler'));
 assert.ok(control(1,'links').querySelector('.cluster-overview'));
 assert.equal(run('JSON.stringify(exportState())'),decisions,'Display controls never approve or edit manual decisions');
 // Click from the extra pane opens the corresponding shared cards, while its
 // own selected locus, source and custom subset remain intact.
 const settings=run('JSON.stringify(comparisons[1].settings())');
 control(1,'links').querySelector('.gene-arrow').dispatchEvent(new w.MouseEvent('click',{bubbles:true}));
 assert.equal(run('JSON.stringify(comparisons[1].settings())'),settings);
 assert.ok(d.querySelector('.gene-focus'),'Gene arrows still expose shared model details');
 // Neighbor review checkboxes intentionally share one decision across viewers.
 change(0,'other','all');change(1,'modelSource','all');change(1,'other','all');
 const a=control(0,'pairTable').querySelector('input[data-pair-key]');
 if(a){a.checked=true;a.dispatchEvent(new w.Event('change'));const b=[...control(1,'pairTable').querySelectorAll('input')].find(n=>n.dataset.pairKey===a.dataset.pairKey);if(b)assert.equal(b.checked,true);assert.equal(run('Object.keys(state.pairs).length'),1)}
 const saved=run('JSON.stringify(comparisons.map(c=>c.settings()))');
 const storage={};for(let i=0;i<w.localStorage.length;i++){const k=w.localStorage.key(i);storage[k]=w.localStorage.getItem(k)}
 const reloaded=open(storage),rw=reloaded.window;
 assert.equal(rw.eval('JSON.stringify(comparisons.map(c=>c.settings()))'),saved,'Pane settings survive a reload independently of the decisions draft');
 assert.equal(rw.eval('JSON.stringify(exportState())'),run('JSON.stringify(exportState())'),'Shared manual decisions survive reload');
 reloaded.window.close();
 control(1,'removeComparison').click();assert.equal(d.querySelectorAll('.comparison-pane').length,2);
 w.dispatchEvent(new w.Event('resize'));assert.ok(control(1,'links').querySelector('.coordinate-ruler'));
 d.getElementById('addComparison').click();assert.equal(d.querySelectorAll('.comparison-pane').length,3);
 assert.equal(run('JSON.stringify(DATA)'),before,'Evidence and model sequences are immutable');
 assert.deepEqual(errors,[]);dom.window.close();
 return {status:'passed',stacked_panes:3,custom_subsets:true,independent_controls:true,shared_decisions:true,reload:true,evidence_unchanged:true};
}
module.exports=checkComparisons;
if(require.main===module)console.log(JSON.stringify(checkComparisons(process.argv[2])));
