const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {webcrypto,createHash}=require('node:crypto'),{JSDOM}=require('jsdom');
const hash=s=>createHash('sha256').update(s).digest('hex');
function data(){
 const proteins={a:'MAAAA',b:'MGGGG'},models=Object.keys(proteins).map(id=>({model_id:id,species:'Synthetic organism',
  original_model:'synthetic_'+id,source:'annotation',protein_aa:5,protein_sha256:hash(proteins[id]),
  contig:'synthetic',start:0,end:18,strand:'+',copy_identity:{display_name:'TEST',basis:'software_control',status:'test'},
  assignment:{best_target:{},facts:['Synthetic software control; no scientific approval.']}}));
 const spec={review_id:'synthetic-reference-ui',set_id:'synthetic-set',protein_sha256:Object.fromEntries(Object.entries(proteins).map(([id,seq])=>[id,hash(seq)])),
  buckets:[{key:'bucket',gene:'TEST',group:'g',genome_id:'synthetic-genome',species:'Synthetic organism',model_ids:['a','b']}]};
 return {groups:[{id:'g',name:'TEST'}],models,finalization:{spec,state:{schema:'stage3_final_choices_v1',review_id:spec.review_id,set_id:spec.set_id,
  choices:{bucket:{primary:'a',reason:'Proposed synthetic representative',reviewed:false}},cases:{}},rows:{primary:[],cases:[]},
  raw_sequences:{protein:proteins,cds:{a:'ATGGCTGCTGCTGCTTAA',b:'ATGGGTGGTGGTGGTTAA'}}},prior_comment_files:{},previous_comments:{},query_reference_files:{},query_protein_count:0};
}
function html(payload=data()){
 const ids=['finalMessage','decisionSearch','decisionCount','decisionTable','specialCaseTable','primaryCount','caseCount','unassignedTargets',
  'finalEditor','finalEditorTitle','finalModelLabel','finalReason','finalEvidence','finalSave','finalCancel','downloadFinal','downloadFinalTop','downloadDefaults','exportFinalJSON','importFinalJSON'];
 const dom=ids.map(id=>id==='finalSave'?'<button id="finalSave">Save choice</button>':id==='finalCancel'?'<button id="finalCancel">Cancel</button>':id==='finalReason'?'<textarea id="'+id+'"></textarea>':id==='decisionSearch'?'<input id="'+id+'">':'<div id="'+id+'"></div>').join('')+'<select id="finalModel"></select>';
 const helpers=`const DATA=${JSON.stringify(payload)},modelMap=new Map(DATA.models.map(m=>[m.model_id,m]));const $=id=>document.getElementById(id);function group(){return 'g'}function pct(x){return x}function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n}function link(href,text,cls){const n=el('a',text,cls);n.href=href;return n}`;
 return '<!doctype html><html lang="ru"><meta charset="utf-8"><title>Проверка кнопки эталона</title><style>body{font:16px/1.6 system-ui;margin:30px;color:#253c41}table{border-collapse:collapse}td,th{padding:14px;border-bottom:1px solid #ddd}button{padding:10px;margin:4px}small,strong,span{display:block}.reference-approval-controls{display:flex;flex-direction:column}#finalMessage{padding:14px;background:#eef6f2}#finalEditor{margin:20px}</style><h1>Проверка интерфейса · шаг 3</h1><p>Синтетические последовательности. Это проверка кнопок, а не реальные модели исследования.</p>'+dom+'<script>'+helpers+'\n'+fs.readFileSync('genecontext/stage3_reference_approval.txt','utf8')+'\n'+fs.readFileSync('genecontext/stage3_finalization.txt','utf8')+'\nfinalChoices();</script></html>';
}
function load(payload=data(),saved){
 const downloads=[];const dom=new JSDOM(html(payload),{url:'http://synthetic-reference.invalid/report.html',runScripts:'dangerously',beforeParse(w){
  w.TextEncoder=TextEncoder;Object.defineProperty(w,'crypto',{value:webcrypto});w.HTMLElement.prototype.scrollIntoView=function(){};
  w.fetch=()=>{throw Error('No prediction submissions or network are allowed in this control')};
  if(saved)w.localStorage.setItem('gcf.final.v1.'+payload.finalization.spec.set_id,saved);
 }});dom.window.__capture=(blob,name)=>downloads.push({blob,name});dom.window.eval('saveBlob=(blob,name)=>window.__capture(blob,name)');
 return {dom,w:dom.window,d:dom.window.document,run:s=>dom.window.eval(s),downloads};
}
const tick=()=>new Promise(r=>setTimeout(r,15));
const readBlob=(w,b)=>new Promise((resolve,reject)=>{const reader=new w.FileReader();reader.onload=()=>resolve(reader.result);reader.onerror=reject;reader.readAsText(b)});
(async()=>{
 const app=load(),before=app.run('JSON.stringify(finalState)');
 assert.equal(app.d.querySelector('[data-action="approve-reference"]'),null,'Old comparison table is not the new STOP3');
 assert.equal(app.d.querySelector('[data-action="send-alphafold"]'),null);
 assert.equal(app.run('JSON.stringify(finalState)'),before);
 app.dom.window.close();console.log(JSON.stringify({status:'passed',legacy_buttons_removed:true,legacy_choices_preserved:true}));
})().catch(e=>{console.error(e);process.exit(1)});
