/* Real browser + compiled launcher + report from a clean wheel installation. */
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const {spawn}=require('node:child_process'),{chromium}=require('playwright');
const root=path.resolve(process.argv[2]),proof=JSON.parse(fs.readFileSync(path.join(root,'verification.json')));
const app=proof.portable,bin=path.join(app,'.review','bin',fs.readdirSync(path.join(app,'.review','bin')).find(n=>n.startsWith('review-')));
const seed=JSON.parse(fs.readFileSync(path.join(app,'review_data.json'))),prefix='gcf-realign-selection-'+seed.review_id;
let child,browser;const errors=[],external=[],submissions=[];
async function start(){
 child=spawn(bin,['--root',app,'--no-browser'],{stdio:['ignore','pipe','pipe']});
 return new Promise((resolve,reject)=>{let out='';const timer=setTimeout(()=>reject(Error('Launcher startup timeout')),15000);child.once('error',reject);child.stdout.on('data',data=>{out+=data;const match=out.match(/Open: (http:\/\/127\.0\.0\.1:[0-9]+\/\?review_token=[a-f0-9]+)/);if(match){clearTimeout(timer);resolve(match[1])}});child.stderr.on('data',data=>{out+=data});child.once('exit',code=>{clearTimeout(timer);reject(Error('Launcher exited: '+code+' '+out))})});
}
async function open(url,first){
 const context=await browser.newContext({viewport:{width:1440,height:1100}});
 if(first)await context.addInitScript(({prefix})=>{
  window.localStorage.setItem('gcf-realign-selection-other-review-saved-alignment-migration-control','old software checkpoint');
  window.localStorage.setItem('gcf-stage2-another-gene','preserve another gene');
 },{prefix});
 const page=await context.newPage();page.on('pageerror',e=>errors.push(e.message));
 page.on('request',r=>{if(r.url().includes('/api/align'))submissions.push(r.url());if(!r.url().startsWith('http://127.0.0.1:')&&!r.url().startsWith('data:'))external.push(r.url())});
 await page.goto(url);await page.waitForSelector('html[data-gcf-review-ready="true"]');
 await page.waitForFunction(()=>window.GCFPortableSession&&window.GCFAlignmentWorkflow);
 return {context,page};
}
async function stop(page){
 await page.getByRole('button',{name:'Finish and clean up',exact:true}).click();
 await page.waitForFunction(()=>document.getElementById('realignStatus').textContent.includes('Local runner closed'));
 await new Promise((resolve,reject)=>{if(child.exitCode!==null)return resolve();const timer=setTimeout(()=>reject(Error('Shutdown timeout')),15000);child.once('exit',()=>{clearTimeout(timer);resolve()})});
}
(async()=>{
 browser=await chromium.launch({headless:true});
 let url=await start(),view=await open(url,true),page=view.page;
 assert.equal(submissions.length,0,'Opening the package must not start a calculation');
 assert.equal(await page.getAttribute('html','data-gcf-checkpoint-storage'),'indexeddb');
 assert.equal(await page.evaluate(key=>window.localStorage.getItem(key),'gcf-realign-selection-other-review-saved-alignment-migration-control'),null);
 assert.equal(await page.evaluate(key=>window.GCFStorage.getItem(key),'gcf-realign-selection-other-review-saved-alignment-migration-control'),'old software checkpoint');
 assert.equal(await page.locator('#syntenyReview').count(),1);
 assert.equal(await page.locator('#initialModelChoices').count(),1);
 assert.equal(await page.evaluate(()=>JSON.parse(document.getElementById('interactiveReviewData').textContent).portable_available),true);
 await page.locator('#alignmentAdvanced').evaluate(e=>e.open=true);
 for(const row of seed.sequences.filter(s=>!s.reference)){const choice=page.getByLabel('Alignment membership '+row.id,{exact:true});await choice.uncheck();await choice.check()}
 await page.getByRole('button',{name:'Save models and build alignment',exact:true}).click();
 await page.waitForSelector('.interactive-aln-row',{timeout:90000});
 await page.waitForFunction(()=>!window.GCFAlignmentWorkflow.current().running);
 assert.equal(submissions.length,1);
 assert.equal(await page.locator('.interactive-aln-label .gcf-alignment-label').count(),2);
 const id=await page.evaluate(()=>window.GCFAlignmentWorkflow.savedAlignment().id);
 await page.locator('#alignmentReviewComment').fill('A portable note tied to this exact alignment');
 await page.evaluate(async({prefix})=>{
  // Capacity control, never read as a scientific alignment by the active group.
  window.GCFStorage.setItem(prefix+'-saved-alignment-software-capacity-control','x'.repeat(18*1024*1024));
  await window.GCFPortableSession.flush();
 },{prefix});
 assert.equal(await page.evaluate(()=>window.localStorage.getItem('gcf-stage2-another-gene')),'preserve another gene');
 await page.evaluate(()=>{window.GCFSyntenyReview.open(false);document.getElementById('syntenyReview').open=false;window.scrollTo(0,0)});
 await page.screenshot({path:path.join(root,'review-desktop.png'),fullPage:true});
 await page.reload();await page.waitForSelector('html[data-gcf-review-ready="true"]');
 assert.equal(await page.evaluate(()=>window.GCFAlignmentWorkflow.savedAlignment().id),id);
 assert.equal(await page.locator('#alignmentReviewComment').inputValue(),'A portable note tied to this exact alignment');
 assert.equal(submissions.length,1,'Reload must reuse the completed result');
 await stop(page);await view.context.close();
 const state=JSON.parse(fs.readFileSync(path.join(app,'.review','session-state.json')));
 assert.equal(state.values[prefix+'-saved-alignment-software-capacity-control'].length,18*1024*1024);
 const previous=url;url=await start();assert.notEqual(url,previous);
 view=await open(url,false);page=view.page;
 assert.equal(await page.evaluate(()=>window.GCFAlignmentWorkflow.savedAlignment().id),id);
 assert.equal(await page.locator('#alignmentReviewComment').inputValue(),'A portable note tied to this exact alignment');
 assert.equal(await page.evaluate(key=>window.GCFStorage.getItem(key).length,prefix+'-saved-alignment-software-capacity-control'),18*1024*1024);
 assert.equal(submissions.length,1,'A fresh browser origin must restore from the package folder');
 await page.evaluate(async key=>{window.GCFStorage.removeItem(key);await window.GCFPortableSession.flush()},prefix+'-saved-alignment-software-capacity-control');
 await page.setViewportSize({width:390,height:844});await page.evaluate(()=>window.scrollTo(0,0));await page.screenshot({path:path.join(root,'review-narrow.png'),fullPage:true});
 await stop(page);await view.context.close();
 assert.deepEqual(errors,[]);assert.deepEqual(external,[]);
 const result={status:'passed',synthetic_only:true,real_mafft:true,indexeddb:true,legacy_migration:true,capacity_mib:18,reload:true,new_origin_restart:true,package_folder_restore:true,no_duplicate_calculations:true,other_gene_unchanged:true,page_errors:errors};
 fs.writeFileSync(path.join(root,'browser-verification.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
})().catch(error=>{console.error(error);console.error('Browser errors:',errors);process.exitCode=1}).finally(async()=>{if(child&&child.exitCode===null)child.kill();await browser?.close()});
