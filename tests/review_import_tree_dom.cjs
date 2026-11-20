const fs=require('node:fs');
const assert=require('node:assert/strict');
const {webcrypto}=require('node:crypto');
const {JSDOM,VirtualConsole}=require('jsdom');

const [page,fixture,profile]=process.argv.slice(2);
const html=fs.readFileSync(page,'utf8');
const example=JSON.parse(JSON.parse(fs.readFileSync(fixture,'utf8')).payload);
const [group]=Object.keys(example.analyses);
const {alignment,tree}=example.analyses[group];
const source=example.source_review_id, review=example.alignment_review_id;
const root=tree.phylogeny.outgroup_sequence_ids[0];
const selection='gcf-realign-selection-'+review;
const oldBasis=tree.phylogeny.root_basis;
const pathname=profile==='hosted'?'/reports/Q0024/report.html':'/report.html?review_token=fixture';
const errors=[];
const requests=[];
const downloads=[];
const wait=()=>new Promise(resolve=>setTimeout(resolve,50));
const equalIDs=(a,b)=>assert.deepEqual([...a].sort(),[...b].sort());
const reviewStorage='gcf-stage2-'+source;
const rootStorage='gcf-outgroup-'+review+'-'+group;
const nativeTree=structuredClone(tree);
nativeTree.phylogeny.mask.trimming={method:'none'};
let rebuiltTree;

function browser(path, storage={}){
 const vc=new VirtualConsole();
 vc.on('jsdomError',e=>{if(!e.message.includes('navigation'))errors.push(e.message)});
 return new JSDOM(html,{url:'http://127.0.0.1:8795'+path,runScripts:'dangerously',virtualConsole:vc,beforeParse(w){
  Object.defineProperty(w,'crypto',{value:webcrypto});w.TextEncoder=TextEncoder;
  w.HTMLElement.prototype.scrollIntoView=function(){};
  Object.entries(storage).forEach(([key,value])=>w.localStorage.setItem(key,value));
  w.URL.createObjectURL=blob=>{
   const reader=new w.FileReader();reader.onload=()=>downloads.push(reader.result);
   reader.readAsText(blob);return 'blob:fixture';
  };
  w.URL.revokeObjectURL=()=>{};w.HTMLAnchorElement.prototype.click=function(){};
  w.fetch=async(url,options)=>{
   const endpoint=new URL(url,w.location.href), action=endpoint.pathname.split('/').pop();
   const reply=value=>({ok:true,json:async()=>JSON.parse(JSON.stringify(value))});
   if(action==='status')return reply({review_id:review,token:'fixture'});
   if(action==='align'){
    const request=JSON.parse(options.body);
    requests.push(request);
    const result=structuredClone(nativeTree);
    result.id='tree-rebuilt-fixture';
    result.phylogeny.outgroup_sequence_ids=request.outgroup_ids;
    result.phylogeny.root_basis=request.root_basis;
    result.phylogeny.mask.trimming={method:request.trim_mode||'none'};
    if(!request.outgroup_ids.length)delete result.phylogeny.rooted_phylogram_svg;
    rebuiltTree=result;
    return reply(profile==='hosted'?{status:'done',ids:request.ids,result}:result);
   }
   if(action==='checkpoint'){
    const id=endpoint.searchParams.get('id');
    return reply(id===alignment.id?alignment:rebuiltTree);
   }
   if(action==='session'||action==='session-state')return reply({values:{}});
   throw new Error('Unexpected request '+endpoint.pathname);
  };
 }});
}

async function upload(w, draft){
 const input=w.document.querySelector('[aria-label="Import alignment model choices"]');
 Object.defineProperty(input,'files',{configurable:true,value:[{text:async()=>JSON.stringify(draft)}]});
 input.dispatchEvent(new w.Event('change',{bubbles:true}));
 await wait();
 return w.document.getElementById('realignStatus').textContent;
}

(async()=>{
 let dom=browser(pathname),w=dom.window;
 await wait();
 assert.notEqual(source,review);
 const initial=w.localStorage.getItem(reviewStorage);
 assert.match(await upload(w,{...example.stage2,review_id:review}),/another STOP 2 review/);
 assert.equal(w.localStorage.getItem(reviewStorage),initial);
 const bad=structuredClone(example.stage2);
 bad.models.unknown={decision:'inspect',role:'target',notes:''};
 assert.match(await upload(w,bad),/lacks 1 selected models/);
 assert.equal(w.localStorage.getItem(reviewStorage),initial);
 assert.match(await upload(w,example.stage2),/Imported 4 CDS models/);
 const draft=JSON.parse(w.localStorage.getItem(reviewStorage));
 assert.deepEqual(draft.models,example.stage2.models);
 assert.deepEqual(draft.regions,example.stage2.regions);
 assert.deepEqual(draft.review_plan,example.stage2.review_plan);
 equalIDs(w.GCFAlignmentWorkflow.current().ids,alignment.ids);
 assert.equal(requests.length,0);
 dom.window.close();

 const stash={
  [reviewStorage]:JSON.stringify(draft),
  [selection]:JSON.stringify(alignment.ids),
  [selection+'-saved-alignment-'+group]:JSON.stringify(alignment),
  [selection+'-input-models-'+alignment.id]:JSON.stringify(alignment.ids),
  [selection+'-result-'+group]:JSON.stringify(nativeTree),
  [rootStorage]:JSON.stringify({ids:[root],basis:oldBasis,trim_al:false}),
 };
 dom=browser(pathname,stash);w=dom.window;
 await wait();
 const doc=w.document;
 const visibleTree=()=>!doc.getElementById('candidateTree').hidden;
 const build=doc.getElementById('buildSelectedTree');
 const pick=doc.querySelector('[aria-label="Outgroup sequences"]');
 const basis=doc.getElementById('treeOutgroupBasis');
 assert.ok(visibleTree());
 const rootRow=[...doc.querySelectorAll('[data-tree-member]')].find(x=>x.dataset.treeMember===root);
 assert.ok(rootRow.checked);
 rootRow.click();
 assert.ok([...pick.options].some(o=>o.value===root&&o.selected));
 assert.ok(!w.GCFAlignmentWorkflow.treeIDs().includes(root));
 assert.ok(visibleTree());
 basis.value='';basis.dispatchEvent(new w.Event('input'));
 assert.ok(build.disabled&&!visibleTree());
 basis.value='Updated software fixture outgroup';basis.dispatchEvent(new w.Event('input'));
 assert.ok(!visibleTree()&&!build.disabled);
 assert.match(doc.getElementById('treeBuildStatus').textContent,/Tree settings changed/);
 doc.getElementById('openSyntenyStop').click();
 assert.equal(w.localStorage.getItem('gcf-synteny-final-'+source),null);
 build.click();
 for(let n=0;n<12&&!visibleTree();n++)await wait();
 assert.ok(visibleTree());
 assert.equal(requests.length,1);
 assert.ok(requests[0].ids.includes(root));
 assert.equal(requests[0].root_basis,basis.value);
 assert.equal(requests[0].alignment_result_id,alignment.id);
 doc.getElementById('openSyntenyStop').click();
 for(let n=0;n<10&&!w.localStorage.getItem('gcf-synteny-final-'+source);n++)await wait();
 assert.equal(JSON.parse(w.localStorage.getItem('gcf-synteny-final-'+source)).tree_result_id,'tree-rebuilt-fixture');
 doc.getElementById('clearOutgroup').click();
 assert.ok(!visibleTree());
 assert.ok(doc.querySelector('[aria-label="Tree rooting view"] option[value="rooted"]').disabled);
 const savedState={};
 for(let i=0;i<w.localStorage.length;i++){
  const key=w.localStorage.key(i);
  savedState[key]=w.localStorage.getItem(key);
 }
 dom.window.close();

 dom=browser(pathname+(pathname.includes('?')?'&':'?')+'view=synteny',savedState);
 w=dom.window;
 const pane=w.document.getElementById('stop4Handoff');
 const roles=[...pane.querySelectorAll('select')];
 assert.equal(roles.length,4);
 for(const role of roles){
  role.value='primary_target';role.dispatchEvent(new w.Event('change',{bubbles:true}));
 }
 pane.querySelector('button').click();
 for(let n=0;n<12&&!downloads.length;n++)await wait();
 assert.equal(downloads.length,1,pane.textContent);
 const payload=JSON.parse(JSON.parse(downloads[0]).payload);
 assert.equal(payload.analyses[group].tree.id,rebuiltTree.id);
 assert.equal(payload.stage3.reviewed,true);
 assert.ok(Object.values(payload.stage3.decisions).every(v=>v.reason===''));
 dom.window.close();
 assert.deepEqual(errors,[]);
 console.log(profile+': import, saved draft, outgroup from an unchecked row, stale-tree guard and STOP 4 export verified');
})().catch(e=>{console.error(e,errors);process.exit(1)});
