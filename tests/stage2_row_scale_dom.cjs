/* Focused linear axes and isolated overlap lanes; exact native CDS coverage. */
const fs=require('fs'),assert=require('node:assert/strict'),crypto=require('crypto');
function check(report,{JSDOM,VirtualConsole}=require('jsdom'),fast=false){
const html=fs.readFileSync(report,'utf8'),start=html.indexOf('const DATA=')+11,end=html.indexOf(',$=id=>',start),data=JSON.parse(html.slice(start,end));
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(html.slice(0,start)+'__exactData'+html.slice(end),{runScripts:'dangerously',pretendToBeVisual:true,url:'http://localhost/qa',virtualConsole:vc,beforeParse(w){w.__exactData=data;w.HTMLElement.prototype.scrollIntoView=function(){}}});
const w=dom.window,d=w.document,run=s=>w.eval(s),draft=run('JSON.stringify(exportState())'),fingerprint=()=>crypto.createHash('sha256').update(JSON.stringify(data)).digest('hex'),before=fingerprint();
let ticks=0,blocks=0,spans=0,lanes=0;const bases=new Map();
function inspect(){
 assert.deepEqual(errors,[]);const svg=d.querySelector('#links .cluster-overview');assert.ok(svg);assert.equal(svg.dataset.scaleMode,'per-row-focused');
 for(const axis of svg.querySelectorAll('.cluster-row-ruler')){
  const row=axis.dataset.row,px=Number(axis.dataset.pxPerBp),center=Number(axis.dataset.rulerCenter),r=data.regions.find(r=>r.region_id===axis.dataset.regionId),mid=(r.anchor_start+r.anchor_end)/2,sign=r.strand==='-'?-1:1;
  const lo=Number(axis.dataset.focusStart),hi=Number(axis.dataset.focusEnd);assert.ok(hi>lo);assert.ok(hi-lo<=Number(axis.dataset.windowBp)+.001);
  for(const tick of axis.querySelectorAll('[data-ruler-bp]')){assert.ok(Math.abs(Number(tick.getAttribute('x1'))-center-Number(tick.dataset.rulerBp)*px)<1e-6);ticks++}
  for(const shape of svg.querySelectorAll(`.gene-arrow[data-row="${row}"] .cluster-model-span`)){
   const a=Number(shape.dataset.start),b=Number(shape.dataset.end),left=center+Math.max(Math.min((a-mid)*sign,(b-mid)*sign),lo)*px,right=center+Math.min(Math.max((a-mid)*sign,(b-mid)*sign),hi)*px;
   assert.ok(Math.abs(Number(shape.getAttribute('x1'))-left)<1e-6);assert.ok(Math.abs(Number(shape.getAttribute('x2'))-right)<1e-6);assert.equal(shape.tagName,'line','Introns must never be filled polygons');spans++;
  }
  for(const shape of svg.querySelectorAll(`.cluster-cds-block[data-row="${row}"]`)){
   const a=Number(shape.dataset.start),b=Number(shape.dataset.end),left=center+Math.max(Math.min((a-mid)*sign,(b-mid)*sign),lo)*px,right=center+Math.min(Math.max((a-mid)*sign,(b-mid)*sign),hi)*px;
   assert.ok(Math.abs(Number(shape.getAttribute('x'))-left)<1e-6);assert.ok(Math.abs(Number(shape.getAttribute('width'))-(right-left))<1e-6);blocks++;
  }
  const groups=[...svg.querySelectorAll(`.gene-arrow[data-row="${row}"]`)];assert.equal(groups.length,svg.querySelectorAll(`.cluster-gene-label[data-row="${row}"]`).length);
  for(let i=0;i<groups.length;i++)for(let j=i+1;j<groups.length;j++)if(groups[i].dataset.lane===groups[j].dataset.lane){assert.ok(Number(groups[i].dataset.layoutRight)<=Number(groups[j].dataset.layoutLeft)+1e-6||Number(groups[j].dataset.layoutRight)<=Number(groups[i].dataset.layoutLeft)+1e-6,'Overlapping genes AND labels must use different lanes');lanes++}
  assert.equal(svg.querySelectorAll(`.cluster-gene-label[data-row="${row}"]`).length,groups.length);
  const key=r.region_id+'|'+d.getElementById('genomeZoom').value+'|'+d.getElementById('links').clientWidth;
  if(d.getElementById('neighborsShown').value==='10'){if(bases.has(key))assert.equal(px,bases.get(key),'Source filters cannot change coordinates');else bases.set(key,px)}
 }
}
const change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
inspect();if(!fast){for(const source of ['annotation','tiberius','all']){change('modelSource',source);inspect()}
for(const z of ['2','0.5','1']){change('genomeZoom',z);inspect()}
change('neighborsShown','5');inspect();change('neighborsShown','10');
for(const width of [390,600,1100]){Object.defineProperty(d.getElementById('links'),'clientWidth',{configurable:true,value:width});run('mainComparison.compare()');inspect();assert.ok(Number(d.querySelector('.cluster-overview').getAttribute('viewBox').split(' ')[2])<=width,'Fit must fit the actual viewer');}
}
assert.equal(fingerprint(),before);assert.equal(run('JSON.stringify(exportState())'),draft);assert.deepEqual(errors,[]);
const result={status:'passed',report,rows:bases.size,ticks,blocks,spans,nonoverlapping_lane_pairs:lanes,scientific_payload_unchanged:true,decisions_unchanged:true,errors};w.close();return result;
}
module.exports=check;if(require.main===module){try{console.log(JSON.stringify(check(process.argv[2],undefined,process.argv.includes('--fast'))))}catch(e){console.error(e.stack);process.exit(1)}}
