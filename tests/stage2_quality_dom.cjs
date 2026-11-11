const fs=require('node:fs'),assert=require('node:assert/strict');
module.exports=function(report,{JSDOM,VirtualConsole}=require('jsdom')){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(report,'utf8'),{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc}),w=dom.window,d=w.document;
 const evidence=w.eval('JSON.stringify(DATA.models)'),review=w.eval('DATA.review_id');
 const change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
 change('region','locus_a327fb5153453624');change('modelSource','all');change('scale','distance');
 assert.ok(d.querySelector('[data-split-prediction="true"]'),'Separate predicted ORFs must be marked on the diagram');
 const recovered=w.eval("!!models.get('Md820f848c1e911c7eb9c6456').reconstruction?.resolves_boundary_conflict");
 assert.ok(d.body.textContent.includes(recovered?'Full CDS reconstructed: 462 aa':'Conflicting gene boundaries: one annotated transcript versus 2 Tiberius ORFs'),'Model state must be visible outside collapsed details');
 assert.ok(d.body.textContent.includes('ORF intact'));assert.ok(!d.body.textContent.includes('Complete CDS'));
 // Focus changes reset vertical pan; zoom preserves it.
 d.getElementById('links').scrollTop=99;change('genomeZoom','2');assert.equal(d.getElementById('links').scrollTop,99);
 change('region','locus_253b72275c40b43e');assert.equal(d.getElementById('links').scrollTop,0);
 assert.ok(w.eval("modelDetails(models.get('Mf532642a0c4c04189f42a18a')).textContent").includes('unclassified translation discrepancy'));
 w.eval(`state.regions=Object.fromEntries(DATA.regions.map(r=>[r.region_id,{decision:r.region_id==='locus_3b4d10c17adaddb3'?'inspect':'exclude',notes:'Isolated software control',saved:true}]));state.models={'M9987e5a0d18a8acbe57af563':{decision:'inspect',role:'target',notes:'Isolated software control',target_override:false}};state.pairs={}`);
 const blocked=w.eval('reviewIssues()');assert.equal(blocked.length,1);assert.ok(blocked[0].message.includes('override'));
 w.eval("state.models.M9987e5a0d18a8acbe57af563.target_override=true");assert.equal(w.eval('reviewIssues().length'),0);
 w.eval(`state.regions=Object.fromEntries(DATA.regions.map(r=>[r.region_id,{decision:r.region_id==='locus_a327fb5153453624'?'inspect':'exclude',notes:'Isolated software control',saved:true}]));state.models={'Md820f848c1e911c7eb9c6456':{decision:'inspect',role:'target',notes:'Isolated software control',target_override:false}}`);
 const annotationBlocked=w.eval('reviewIssues()');assert.equal(annotationBlocked.length,recovered?0:1);if(!recovered)assert.ok(annotationBlocked[0].message.includes('Conflicting gene boundaries'));
 w.eval("state.models.Md820f848c1e911c7eb9c6456.target_override=true");assert.equal(w.eval('reviewIssues().length'),0);
 assert.equal(w.eval('JSON.stringify(DATA.models)'),evidence);assert.equal(w.eval('DATA.review_id'),review);assert.deepEqual(errors,[]);
 dom.window.close();return {status:'passed',split_visible:true,short_model_guard:true,annotation_conflict_state:recovered?'reference_backed_full_model':'requires_override',explicit_override:true,annotation_flags:true,scientific_models_unchanged:true};
};
