// The package folder is the portable source of saved review state, including
// IndexedDB-backed checkpoints. Restore it before any report scripts execute.
window.GCFPortableReady=(async()=>{
 'use strict';
 const config=JSON.parse(document.getElementById('gcfPortableState').textContent);
 const state=config.state,token=new URLSearchParams(location.search).get('review_token');
 if(!token)return;
 await window.GCFStorageReady;
 const store=window.GCFStorage||window.localStorage;
 const prefixes=['gcf-realign-selection-'+config.review_id,'gcf-review-workspace-2-'+config.review_id,'gcf-review-workspace-3-'+config.review_id,'gcf-one-step-undo-'+config.review_id,'gcf-outgroup-'+config.review_id,'gcf-alignment-first-notes-'+config.review_id,...(config.source_review_id?['gcf-stage2-'+config.source_review_id,'gcf-stage2-views-'+config.source_review_id,'gcf-synteny-final-'+config.source_review_id,'gcf-synteny-reference-v1-'+config.source_review_id,'gcf-stop4-final-'+config.source_review_id]:[])];
 const allowed=k=>prefixes.some(p=>k===p||k.startsWith(p+'-'));
 const keys=()=>store.keys?store.keys():Array.from({length:store.length},(_,i)=>store.key(i));
 // A new port/origin must restore the same state, without touching another gene.
 for(const key of keys())if(allowed(key))store.removeItem(key);
 for(const[key,value]of Object.entries(state.values||{}))if(allowed(key))store.setItem(key,value);
 await store.flush?.();
 let dirty=false,saving=false,timer,lastError='',inFlight;
 function banner(){const el=document.getElementById('portableSaveStatus');if(el)el.textContent=lastError?'Not saved to the package folder: '+lastError:saving?'Saving to the package folder…':dirty?'Unsaved changes':'Choices and results saved to the package folder'}
 async function saveState(){
  if(!dirty)return;saving=true;dirty=false;banner();
  try{
   await store.flush?.();
   const values=Object.fromEntries(keys().filter(allowed).map(k=>[k,store.getItem(k)]));
   const response=await fetch('/api/state',{method:'POST',headers:{'Content-Type':'application/json','X-Review-Token':token},body:JSON.stringify({schema:state.schema,dataset_sha256:state.dataset_sha256,values})});
   if(!response.ok)throw Error(await response.text());lastError='';
  }catch(error){dirty=true;lastError=error.message}
  finally{saving=false;banner()}
 }
 async function flush(){
  clearTimeout(timer);
  while(dirty||inFlight){if(!inFlight)inFlight=saveState().finally(()=>{inFlight=null});await inFlight;if(lastError)throw Error(lastError)}
  await store.flush?.();
 }
 function changed(key){if(!allowed(String(key)))return;dirty=true;lastError='';clearTimeout(timer);timer=setTimeout(()=>flush().catch(()=>{}),200);banner()}
 if(store.subscribe)store.subscribe(changed);
 else{
  const set=Storage.prototype.setItem,remove=Storage.prototype.removeItem;
  Storage.prototype.setItem=function(k,v){set.call(this,k,v);if(this===store)changed(k)};
  Storage.prototype.removeItem=function(k){remove.call(this,k);if(this===store)changed(k)};
 }
 window.GCFPortableSession={flush,config};
 function showBanner(){if(document.getElementById('portableSaveStatus'))return;const el=document.createElement('p');el.id='portableSaveStatus';el.style='padding:8px 12px;background:#edf5f1;color:#234b42';document.querySelector('main')?.prepend(el);banner()}
 if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',showBanner,{once:true});else showBanner();
 window.addEventListener('beforeunload',event=>{if(dirty||saving){event.preventDefault();event.returnValue='Wait until changes are saved to the package folder.'}});
})();
