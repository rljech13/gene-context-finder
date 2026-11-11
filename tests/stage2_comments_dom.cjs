const fs=require('node:fs'),assert=require('node:assert/strict');
module.exports=function(report,{JSDOM,VirtualConsole}=require('jsdom')){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(report,'utf8'),{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc}),w=dom.window,d=w.document;
 const before=w.eval('JSON.stringify(DATA.models)'),change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
 change('region','locus_a327fb5153453624');change('modelSource','all');change('scale','distance');
 const fragments=w.eval('DATA.filtered_predictions');assert.equal(fragments.length,1);const f=fragments[0];assert.equal(f.cds_length,132);assert.equal(f.target_cds_overlap,112);assert.equal(f.reviewable,false);
 assert.equal(w.eval('models.has(DATA.filtered_predictions[0].model_id)'),false);assert.ok(d.querySelector('.cluster-model-span[data-filtered="true"]'),'Filtered target CDS must be visible');
 assert.equal(d.querySelectorAll('[data-link-kind="miniprot"]').length,0);
 const notes=[...d.querySelectorAll('.automatic-comment')].filter(n=>n.textContent.includes('214.5 kb'));assert.ok(notes.length);for(const n of notes)assert.ok(n.closest('details')&&!n.closest('details').open,'Comments must be collapsed by default');
 const shapes=[...d.querySelectorAll('.cluster-model-span')];assert.ok(shapes.some(s=>s.dataset.sourceBadge==='A'));assert.ok(shapes.some(s=>s.dataset.sourceBadge==='T'));
 change('modelSource','annotation');assert.equal(d.querySelectorAll('.cluster-model-span[data-filtered="true"]').length,0);
 change('modelSource','tiberius');assert.ok(d.querySelector('.cluster-model-span[data-filtered="true"]'));
 w.eval("state.regions['locus_a327fb5153453624']={notes:'Reviewer text kept',decision:'',saved:false}");const exported=w.eval('exportState()');assert.equal(exported.regions['locus_a327fb5153453624'].notes,'Reviewer text kept');assert.ok(exported.automatic_comments['locus_a327fb5153453624'].length>=2);
 assert.equal(w.eval('JSON.stringify(DATA.models)'),before);assert.equal(w.eval('Object.keys(state.models).length'),0);assert.deepEqual(errors,[]);dom.window.close();
 return {status:'passed',annotation_and_tiberius:true,filtered_132_nt_visible:true,filtered_not_selectable:true,automatic_comments_collapsed:true,manual_text_preserved:true,no_approvals:true};
};
