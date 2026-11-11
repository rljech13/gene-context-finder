/* NODE_PATH=/tmp/gcf_report_ui_test/node_modules node tests/stage2_zoom_dom.cjs REPORT */
const fs=require('node:fs'),assert=require('node:assert/strict');
function checkZoom(report,{JSDOM,VirtualConsole}=require('jsdom')){
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(report,'utf8'),{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc});
const w=dom.window,d=w.document,before=JSON.stringify(w.eval('DATA'));
const change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
const results=[];
for(const view of ['distance','models']){
 change('scale',view);change('genomeZoom','1');
 const svg=()=>d.querySelector('#links svg'),scale=()=>Number(svg().dataset.pxPerBp);
 const base=scale(),height=()=>Number(svg().getAttribute('viewBox').split(' ')[3]),baseHeight=height();
 const pane=d.getElementById('links'),focus=()=>(pane.scrollLeft+(pane.clientWidth+250)/2-Number(svg().dataset.rulerCenter))/scale();
 pane.scrollLeft+=120;pane.scrollTop=37;const initialFocus=focus();
 for(const zoom of [2,4,1]){
  change('genomeZoom',String(zoom));
  assert.ok(Math.abs(scale()/base-zoom)<1e-10,`${view}: ${zoom}x must multiply genomic distance by ${zoom}, including wide label-packed canvases`);
  if(view==='distance')assert.equal(height(),baseHeight,'Zoom must not stretch the fixed-height locus tracks');
  assert.ok(Math.abs(focus()-initialFocus)<1e-7,'Zoom must preserve the genomic viewport focus');
  assert.equal(pane.scrollTop,37,'Zoom must preserve vertical review position');
  const marker=view==='distance'?'.cluster-model-span':'.cds-block';
  for(const shape of svg().querySelectorAll(marker)){
   const bp=Number(shape.dataset.end)-Number(shape.dataset.start);
   const xs=shape.getAttribute('points')?.split(/\s+/).map(p=>Number(p.split(',')[0]));
   const width=xs?Math.max(...xs)-Math.min(...xs):Number(shape.getAttribute('width'));
   assert.ok(Math.abs(width-bp*scale())<1e-7,'Zoom must retain native bp lengths');
  }
  results.push({view,zoom,pxPerBp:scale()});
 }
}
// Source filters, focal loci, neighbor counts and comparison subsets must
// not change the meaning of one bp or redraw the reviewed window endpoints.
const scaleChecks=[];
change('scale','distance');change('genomeZoom','1');
const commonScale=Number(d.querySelector('#links svg').dataset.pxPerBp);
assert.equal(parseFloat(d.querySelector('#links svg').style.width),1260,'Fit must use the whole panel without a global horizontal overflow (JSDOM fallback width)');
const checkScale=()=>{
 const svg=d.querySelector('#links svg'),px=Number(svg.dataset.pxPerBp);
 assert.equal(px,commonScale,'One bp must have the same width in every source/locus/comparison view');
 assert.equal(parseFloat(svg.style.width),Number(svg.getAttribute('viewBox').split(' ')[2]),'The browser must not fit-rescale the genomic SVG');
 const rows=w.eval('comparisonRows()');
 for(const [i,r]of rows.entries()){
  const a=svg.querySelector(`[data-window-edge="left"][data-row="${i}"]`),b=svg.querySelector(`[data-window-edge="right"][data-row="${i}"]`);
  assert.ok(Math.abs(Number(b.getAttribute('x1'))-Number(a.getAttribute('x1'))-(r.window_end-r.window_start)*px)<1e-7,'Window width must equal its saved native bp span');
  assert.equal(Number(a.dataset.coordinate),r.strand==='-'?r.window_end:r.window_start);
  assert.equal(Number(b.dataset.coordinate),r.strand==='-'?r.window_start:r.window_end);
  if(d.getElementById('scale').value==='distance'){
   const local=svg.querySelector(`.cluster-row-ruler[data-row="${i}"]`);assert.ok(local,'Every cluster row needs its own aligned kb ruler');
   const ticks=[...local.querySelectorAll('[data-ruler-bp]')];for(let j=1;j<ticks.length;j++){const delta=Number(ticks[j].dataset.rulerBp)-Number(ticks[j-1].dataset.rulerBp);assert.ok(Math.abs(Number(ticks[j].getAttribute('x1'))-Number(ticks[j-1].getAttribute('x1'))-delta*px)<1e-8)}
   assert.equal(svg.querySelectorAll('.target-span-label').length,0,'Original NANOG captions should not gain model-span summary rows');
  }
  for(const label of svg.querySelectorAll(`.cluster-gene-label[data-row="${i}"]`)){
   assert.ok(Number(label.dataset.left)>=0,'Readable labels must remain on the SVG canvas');
   assert.ok(Number(label.dataset.right)<=Number(svg.getAttribute('viewBox').split(' ')[2])+1e-7,'Names may move outside the DNA span, but remain on the canvas');
  }
 }
 scaleChecks.push({view:d.getElementById('scale').value,source:d.getElementById('modelSource').value,pxPerBp:px});
};
for(const view of ['distance','models']){
 change('scale',view);
 for(const source of ['all','annotation','tiberius']){change('modelSource',source);checkScale()}
}
change('scale','distance');change('modelSource','all');
for(const r of w.eval('DATA.regions')){change('region',r.region_id);checkScale()}
Object.defineProperty(d.getElementById('links'),'clientWidth',{configurable:true,value:760});w.dispatchEvent(new w.Event('resize'));assert.equal(parseFloat(d.querySelector('#links svg').style.width),760,'Fit must follow the inner panel width after resizing');Object.defineProperty(d.getElementById('links'),'clientWidth',{configurable:true,value:0});w.dispatchEvent(new w.Event('resize'));checkScale();
change('neighborsShown','5');checkScale();change('neighborsShown','10');checkScale();
const other=[...d.getElementById('other').options].find(o=>o.value!=='all');
if(other){change('other',other.value);checkScale();change('other','all')}
assert.equal(JSON.stringify(w.eval('DATA')),before,'Display zoom must not modify scientific evidence');
assert.equal(w.eval('Object.keys(state.models).length'),0,'Display zoom must not approve models');
assert.deepEqual(errors,[]);
dom.window.close();
return {status:'passed',results,scaleChecks};
}
module.exports=checkZoom;
if(require.main===module)console.log(JSON.stringify(checkZoom(process.argv[2])));
