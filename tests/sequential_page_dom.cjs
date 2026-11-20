/* Exercise the complete current page, including file viewing with no browser storage. */
const fs=require('node:fs'),assert=require('node:assert/strict'),{webcrypto}=require('node:crypto'),{JSDOM,VirtualConsole}=require('jsdom');
const page=fs.readFileSync(process.argv[2],'utf8');
(async()=>{for(const url of ['file:///control/report.html','https://control.invalid/report.html','https://control.invalid/report.html?view=synteny']){
const errors=[],requests=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));const dom=new JSDOM(page,{url,runScripts:'dangerously',pretendToBeVisual:true,virtualConsole:vc,beforeParse(w){w.TextEncoder=TextEncoder;Object.defineProperty(w,'crypto',{value:webcrypto});w.HTMLElement.prototype.scrollIntoView=function(){};w.fetch=async u=>{requests.push(String(u));throw Error('No runner configured in software control')}}});await new Promise(r=>setTimeout(r,40));const d=dom.window.document;assert.deepEqual(errors,[]);
if(url.includes('synteny')){assert.ok(d.documentElement.classList.contains('gcf-synteny-stop'));assert.ok(d.querySelector('#syntenyReferenceApproval'));assert.ok([...d.querySelectorAll('[data-action="send-alphafold"]')].every(b=>b.disabled));}
else{assert.ok(d.querySelector('#initialModelChoices'));assert.ok(d.querySelector('#treeBuildStep'));assert.ok(d.querySelector('#finalContextHandoff'));assert.ok(d.querySelector('a[href="review_app.zip"]').hidden);}
if(url.startsWith('file:'))assert.equal(requests.length,0);assert.equal(requests.filter(u=>u.includes('/api/align')).length,0);dom.window.close();}
console.log(JSON.stringify({status:'passed',complete_page:true,file_without_storage:true,new_stop3:true,no_calculation_on_open:true}))})().catch(e=>{console.error(e.stack);process.exit(1)});
