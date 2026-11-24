const fs=require('node:fs'),assert=require('node:assert/strict');
function check(report,{JSDOM,VirtualConsole}=require('jsdom')){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(report,'utf8'),{runScripts:'dangerously',url:'http://localhost/inherited-comments/',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}});
 const w=dom.window,d=w.document,run=s=>w.eval(s),before=run('JSON.stringify(DATA)'),draft=run('JSON.stringify(exportState())');
 const r=run('DATA.regions[0]');assert.ok(r.notes.includes('KLF2 - has some similarity'),'Use the actual KLF4 selection note as the fixture');
 assert.equal(d.querySelector('#stage1Comment .stop1-comment-body').textContent,r.notes);
 assert.equal(d.querySelector('#stage1Comment textarea'),null,'Inherited notes are read-only');
 const row=d.querySelector(`tr[data-region-id="${r.region_id}"][data-primary="true"]`);
 assert.equal(row.querySelector('.stop1-comment-body').textContent,r.notes);
 assert.equal(row.querySelector('textarea').value,'','An inherited comment is not a STOP-2 review decision');
 d.getElementById('addComparison').click();
 const other=run('DATA.regions.find(r=>r.notes&&r.notes!==DATA.regions[0].notes)');assert.ok(other);
 const sel=d.getElementById('comparison2-region');sel.value=other.region_id;sel.dispatchEvent(new w.Event('change'));
 assert.equal(d.querySelector('#comparison2-stage1Comment .stop1-comment-body').textContent,other.notes);
 assert.equal(d.querySelector('#stage1Comment .stop1-comment-body').textContent,r.notes,'Each viewer follows its own focal locus');
 assert.equal(run('JSON.stringify(exportState())'),draft,'Viewing notes does not approve models or replace the new review comment');
 // Plain text, line breaks and the absence of a recorded comment remain clear.
 w.__note={region_id:'fixture',notes:'First line\n<b>plain text</b>',reason:'other'};
 assert.equal(run('stop1Comment(__note).querySelector("p").textContent'),'First line\n<b>plain text</b>');
 assert.equal(run('stop1Comment(__note).querySelector("b")'),null);
 assert.match(run('stop1Comment({region_id:"empty"},true).textContent'),/No comment recorded/);
 assert.equal(run('stop1Comment({region_id:"empty"}).hidden'),true);
 assert.equal(run('JSON.stringify(DATA)'),before);assert.deepEqual(errors,[]);dom.window.close();
 return {status:'passed',original_comment_exact:true,visible_in_viewer_and_decision_table:true,independent_viewers:true,read_only:true,no_autoapproval:true,plain_text_safe:true};
}
module.exports=check;if(require.main===module)console.log(JSON.stringify(check(process.argv[2])));
