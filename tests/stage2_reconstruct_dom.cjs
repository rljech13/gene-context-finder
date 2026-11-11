const fs=require('node:fs'),assert=require('node:assert/strict');
module.exports=function(report,{JSDOM,VirtualConsole}=require('jsdom')){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(report,'utf8'),{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc}),w=dom.window,d=w.document;
 const change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
 change('region','locus_a327fb5153453624');change('scale','distance');
 const model=w.eval("models.get('Md820f848c1e911c7eb9c6456')"),rec=model.reconstruction;
 assert.equal(rec.cds_length,1389);assert.equal(rec.protein_length,462);assert.equal(rec.transcript_length,2089);assert.equal(rec.sequence_identity,1);assert.equal(rec.matched_bases,2089);assert.equal(rec.cds_blocks.length,5);
 assert.equal(w.eval("modelConstraint(models.get('Md820f848c1e911c7eb9c6456')).target_warnings.length"),0);
 assert.ok(w.eval("modelConstraint(models.get('M314cf2539d20efa20bce69f7')).target_warnings").some(x=>x.includes('Fragmentary prediction')));
 const label=d.querySelector('.cluster-gene-label[data-row="0"]');assert.equal(label.dataset.sourceBadge,'A/T');assert.ok(!label.textContent.includes('full model'));assert.ok(!label.querySelector('text').textContent.includes('family candidate'));assert.ok(label.querySelector('title').textContent.includes('family candidate'));assert.ok(!label.textContent.includes('Unassigned'));
 assert.ok(!d.getElementById('regionInfo').textContent.includes('transcript:'),'Reconstruction metrics must not crowd the diagram header');
 assert.equal(d.getElementById('reconstructionFiles').hidden,false);assert.equal(d.querySelectorAll('#reconstructionFiles a').length,3);
 const menu=d.querySelector('[aria-label="Main model locus_a327fb5153453624"]');assert.ok(menu.selectedOptions[0].textContent.startsWith('Annotation + transcript'));
 assert.equal(w.eval('Object.keys(state.models).length'),0);assert.deepEqual(errors,[]);dom.window.close();
 return {status:'passed',full_cds_length:1389,full_protein_length:462,coding_exons:5,transcript_matched_bases:2089,tiberius_fragments_retained:true,downloads:3,no_automatic_approval:true};
};
