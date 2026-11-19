// Software control: standalone embedded export, defaults vs edited choices.
import fs from 'node:fs';
import {Blob} from 'node:buffer';
import {TextEncoder} from 'node:util';
import {URL} from 'node:url';
import {setTimeout} from 'node:timers';
import vm from 'node:vm';
import assert from 'node:assert/strict';
export async function checkExport(template, archiveDir){
 const models=['a','b'].map(model_id=>({model_id,species:'Species A',genome_id:'A',original_model:model_id,source:'annotation',contig:'chr1',start:0,end:12,strand:'+',protein_aa:3,copy_identity:{display_name:'GENE',status:'source',basis:'annotation'},assignment:{facts:['Test evidence'],best_target:{}}}));
 const state={schema:'stage3_final_choices_v1',review_id:'software',set_id:'control',choices:{bucket:{primary:'a',reason:'Default proposal',reviewed:false}},cases:{}};
 const data={groups:[{name:'GENE'}],models,finalization:{spec:{review_id:'software',set_id:'control',buckets:[{key:'bucket',gene:'GENE',group:'g',model_ids:['a','b']}]},state,rows:{primary:[],cases:[]},raw_sequences:{protein:{a:'MAA',b:'MCC'},cds:{a:'ATGGCTGCTTAA',b:'ATGTGTTGTTAA'}}},prior_comment_files:{'comments/stop1.json':'{"notes":"literal"}\r\n'},previous_comments:{a:{stop1_notes:'Стоп 1\nКомментарий\t"цитата"'},b:{stop2_model_notes:'Копия B'}},comment_history:{status:'software'}};
 data.query_reference_files={'query_references.faa':'>reference_q1\nMQQ\n>reference_q2\nMQQ\n','query_references.tsv':'reference_id\taccession\nreference_q1\tNP_test.1\nreference_q2\tNP_test.2\n'};data.query_protein_count=2;
 const elements=new Map(),store=new Map(),saved=[];
 const ctx=vm.createContext({DATA:data,modelMap:new Map(models.map(m=>[m.model_id,m])),Blob,TextEncoder,setTimeout,URL,localStorage:{getItem:k=>store.get(k),setItem:(k,v)=>store.set(k,v)},$:id=>{if(!elements.has(id))elements.set(id,{listeners:{},addEventListener(event,fn){this.listeners[event]=fn},textContent:''});return elements.get(id)},capture:(blob,name)=>saved.push({blob,name})});
 vm.runInContext(fs.readFileSync(template,'utf8'),ctx);
 vm.runInContext('saveBlob=(blob,name)=>capture(blob,name)',ctx);
 const original=JSON.stringify(data);
 // No row saves: current action must include the full default set.
 elements.get('downloadFinalTop').listeners.click();
 assert.equal(saved.length,1);assert.ok(saved[0].name.endsWith('.chosen_models_with_comments.zip'));
 vm.runInContext("finalState.choices.bucket={primary:'b',reason:'Reviewed alternative\\nwith comment',reviewed:true};persistFinal()",ctx);
 const edited=vm.runInContext('JSON.stringify(finalState)',ctx),stored=[...store.entries()];
 elements.get('downloadFinal').listeners.click();elements.get('downloadDefaults').listeners.click();
 assert.equal(saved.length,3);assert.equal(JSON.stringify(data),original);assert.equal(vm.runInContext('JSON.stringify(finalState)',ctx),edited);assert.deepEqual([...store.entries()],stored);
 const defaults=JSON.parse(vm.runInContext('finalFiles(DATA.finalization.state)["stage3.final.decisions.json"]',ctx));assert.deepEqual(defaults,state);
 const current=JSON.parse(vm.runInContext('finalFiles()["stage3.final.decisions.json"]',ctx));assert.equal(current.choices.bucket.primary,'b');assert.equal(current.choices.bucket.reviewed,true);
 const files=vm.runInContext('finalFiles()',ctx);assert.ok(files['models_with_comments.tsv'].includes('Стоп 1\nКомментарий\t""цитата""'));assert.ok(files['models_with_comments.tsv'].includes('stop2_model_notes'));assert.equal(files['comments/stop1.json'],data.prior_comment_files['comments/stop1.json']);
 assert.equal(files['query_references.faa'],data.query_reference_files['query_references.faa']);assert.equal(files['all_selected_with_queries.proteins.faa'].split('>').length-1,4);assert.equal(files['primary_targets_with_queries.proteins.faa'].split('>').length-1,3);assert.ok(files['primary_targets_with_queries.proteins.faa'].startsWith('>b\n'));assert.equal(defaults.choices.bucket.primary,'a');
 for(const name of ['primary_targets','special_cases','all_selected'])for(const suffix of ['proteins.faa','cds.fna'])assert.ok(name+'.'+suffix in files);
 if(archiveDir){fs.mkdirSync(archiveDir,{recursive:true});for(let i=0;i<saved.length;i++)fs.writeFileSync(archiveDir+'/'+i+'.zip',new Uint8Array(await saved[i].blob.arrayBuffer()))}
 return {passed:true,default_without_row_saves:true,current_and_original_defaults_separate:true,literal_comments_preserved:true,all_native_sequences:true,all_queries_in_both_exports:true,distinct_identical_queries_preserved:true,review_state_unchanged_by_export:true,archives:saved.length};
}
