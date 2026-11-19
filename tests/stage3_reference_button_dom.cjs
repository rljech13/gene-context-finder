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
 let app=load(),{w,d,run}=app;const button=action=>d.querySelector('[data-action="'+action+'"]');
 assert.equal(button('approve-reference').disabled,true);assert.equal(button('send-alphafold').disabled,true);
 run("editPrimary('a')");d.getElementById('finalReason').value='Reviewed synthetic final model';d.getElementById('finalSave').click();
 assert.equal(button('approve-reference').disabled,false);assert.equal(button('send-alphafold').disabled,true,'Saving a reviewed choice does not approve a reference');
 button('approve-reference').click();await tick();
 assert.equal(button('send-alphafold').disabled,false);assert.equal(run('finalState.choices.bucket.reference_approval.protein_sha256'),hash('MAAAA'));
 const approved=run('JSON.stringify(finalState)');assert.ok(JSON.parse(run('localStorage.getItem(finalKey)')).choices.bucket.reference_approval);
 button('send-alphafold').click();await tick();assert.equal(app.downloads.length,1);
 const job=JSON.parse(await readBlob(w,app.downloads[0].blob));assert.equal(job.length,1);assert.equal(job[0].sequences[0].proteinChain.sequence,'MAAAA');assert.deepEqual(job[0].modelSeeds,['42']);
 assert.match(d.getElementById('finalMessage').textContent,/Задание скачано/);assert.equal(d.querySelector('#finalMessage a').href,'https://alphafoldserver.com/');
 assert.ok(JSON.parse(run("finalFiles()['stage3.final.decisions.json']")).choices.bucket.reference_approval);
 app.dom.window.close();app=load(data(),approved);({w,d,run}=app);assert.equal(button('send-alphafold').disabled,false,'Explicit approval survives reload');
 button('revoke-reference').click();assert.equal(button('send-alphafold').disabled,true);assert.equal(run('finalState.reference_history.at(-1).action'),'revoked');
 button('approve-reference').click();await tick();run("editPrimary('b')");d.getElementById('finalReason').value='Changed synthetic final model';d.getElementById('finalSave').click();
 assert.equal(button('send-alphafold').disabled,true);assert.equal(run('finalState.choices.bucket.reference_approval'),undefined);assert.equal(run('finalState.reference_history.at(-1).action'),'invalidated_by_model_edit');
 app.dom.window.close();const changed=data();changed.finalization.raw_sequences.protein.a='MAAAC';changed.models[0].protein_sha256=hash('MAAAC');changed.finalization.spec.protein_sha256.a=hash('MAAAC');app=load(changed,approved);({w,d,run}=app);
 assert.equal(button('send-alphafold').disabled,true,'Changed sequence invalidates old approval');assert.equal(run('finalState.choices.bucket.primary'),'a');assert.equal(run('finalState.choices.bucket.reviewed'),true);assert.equal(run('finalState.reference_history.at(-1).action'),'invalidated_on_load');
 app.dom.window.close();
 if(process.argv[2]){fs.mkdirSync(process.argv[2],{recursive:true});fs.writeFileSync(path.join(process.argv[2],'report.html'),html());fs.writeFileSync(path.join(process.argv[2],'approved.synthetic.json'),approved)}
 console.log(JSON.stringify({status:'passed',legacy_review_does_not_approve:true,explicit_button_required:true,exact_sequence_bound:true,reload:true,revoke:true,changed_model_blocks:true,changed_sequence_blocks:true,one_job_export:true,no_network_or_predictions:true}));
})().catch(e=>{console.error(e);process.exit(1)});
