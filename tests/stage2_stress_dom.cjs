/* Bounded synthetic UI controls. No production evidence/decisions are loaded. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {JSDOM,VirtualConsole}=require('jsdom');
const root=path.resolve(__dirname,'../genecontext'),template=fs.readFileSync(path.join(root,'stage2_review.txt'),'utf8').replace('__NEIGHBORHOOD__',fs.readFileSync(path.join(root,'stage2_neighborhood.txt'),'utf8')).replace('__CLUSTER__',fs.readFileSync(path.join(root,'stage2_cluster.txt'),'utf8')).replace('__REVIEW_QUEUE__',fs.readFileSync(path.join(root,'stage2_review_queue.txt'),'utf8'));
function fixture(count){
 const regions=[0,1,2].map(i=>({region_id:'r'+i,genome_id:i===2?'Species_a_GCF_1.1':i===0?'Species_a_GCF_1.1':'Species_b_GCF_2.1',contig:'chr',strand:i===0?'-':'+',window_start:0,window_end:1000000,anchor_start:500000,anchor_end:501000}));
 const models=[];
 for(const r of regions)for(let k=-count;k<=count;k++)for(const source of ['annotation','tiberius']){
  const g=(r.genome_id===regions[0].genome_id?'a':'b')+'_'+k;
  const start=k===0?500010:500000+k*10000,end=start+500;
  models.push({model_id:r.region_id+'_'+k+'_'+source,region_id:r.region_id,genome_id:r.genome_id,gene_group:g,original_model:'tx_'+k,original_gene:g,source,annotation_name:k===0?'TEST':'GENE_'+k,contig:'chr',strand:k%2?'-':'+',start,end,cds_blocks:[[start,end,0]],stop_blocks:[],overlaps_anchor:k===0,window_edge:false,coding_status:'complete_predicted',issues:[],protein_length:100,identification:{proposed_name:k===0?'TEST':'GENE_'+k,near_best:[],status:'homology_assignment'},target_support:[]});
 }
 const neighbor_pairs=[];
 for(let k=-count;k<=count;k++)neighbor_pairs.push({gene_a:'a_'+k,gene_b:'b_'+k,model_a:'r0_'+k+'_annotation',model_b:'r1_'+k+'_annotation',kind:k===0?'target_similarity':'neighbor_proposal',status:'unique_reciprocal_proposal',identity:.8,coverage_query:1,coverage_subject:1,bitscore:100,evalue:1e-20});
 const source_neighbor_pairs={annotation:neighbor_pairs,tiberius:neighbor_pairs.map(p=>({...p,model_a:p.model_a.replace(/annotation$/,'tiberius'),model_b:p.model_b.replace(/annotation$/,'tiberius')}))};
 return {review_id:'isolated-ui-control',source_review_id:'fixture',status:'awaiting_manual_review',regions,models,neighbor_pairs,source_neighbor_pairs,errors:[],query_members:{q:[{query:'control'}]}};
}
function open(data,suffix=''){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(template.replace('__DATA__',JSON.stringify(data)),{runScripts:'dangerously',url:'http://localhost/isolated-control'+suffix,virtualConsole:vc});
 const w=dom.window,d=w.document,change=(id,value)=>{d.getElementById(id).value=value;d.getElementById(id).dispatchEvent(new w.Event('change'))};
 return {dom,w,d,errors,change,run:s=>w.eval(s)};
}
(async()=>{
 const start=performance.now(),data=fixture(22),dense=open(data),{d,w,run,change}=dense,before=JSON.stringify(run('DATA'));
 await new Promise(r=>setImmediate(r));assert.deepEqual(dense.errors,[]);
 assert.equal(d.querySelectorAll('#models details.gene').length,45);
 let arrows=d.querySelectorAll('#links .gene-arrow[data-row="0"]');assert.equal(arrows.length,21);
 assert.equal(d.querySelector('#links .gene-count[data-row="0"]').dataset.total,'45');
 assert.equal(d.querySelector('#links .gene-count[data-row="0"]').dataset.shown,'21');
 assert.equal(run('neighborhood(current()).hiddenLeft'),12);assert.equal(run('neighborhood(current()).hiddenRight'),12);
 assert.match(d.querySelector('#links [data-window-row="0"]').textContent,/12 hidden/);
 change('neighborsShown','5');assert.equal(d.querySelectorAll('#links .gene-arrow[data-row="0"]').length,11);
 assert.equal(d.querySelectorAll('#models details.gene').length,45,'Display limit must not remove evidence or model choices');
 change('neighborsShown','10');change('other','all');
 assert.equal(Number(d.querySelector('.cluster-overview').getAttribute('viewBox').split(' ')[3]),112*3+68);
 assert.equal(Number(d.querySelector('.cluster-overview').dataset.rowHeight),112);
 for(const label of d.querySelectorAll('.cluster-gene-label')){assert.ok(Number(label.dataset.left)>=249);assert.ok(Number(label.dataset.right)<=Number(d.querySelector('.cluster-overview').getAttribute('viewBox').split(' ')[2])-49)}
 assert.equal(new Set([...d.querySelectorAll('#links .gene-arrow')].map(n=>n.dataset.row)).size,3);
 for(const mode of ['order','models']){
  change('scale',mode);
  if(mode==='models'){
   for(const rect of d.querySelectorAll('#links .cds-block')){
    assert.ok(Number(rect.getAttribute('x'))>=90&&Number(rect.getAttribute('x'))+Number(rect.getAttribute('width'))<=1170);
    assert.ok(Number(rect.getAttribute('width'))>0);
   }
   const labels=[...d.querySelectorAll('#links .ruler-gene-label')];
   for(let i=0;i<labels.length;i++)for(let j=i+1;j<labels.length;j++){
    const a=labels[i].dataset,b=labels[j].dataset;
    assert.equal(Number(a.left)<Number(b.right)&&Number(b.left)<Number(a.right)&&Number(a.top)<Number(b.bottom)&&Number(b.top)<Number(a.bottom),false);
   }
  }

  for(const arrow of d.querySelectorAll('#links polygon.source-arrow'))for(const point of arrow.getAttribute('points').split(/\s+/)){
   const [x,y]=point.split(',').map(Number);assert.ok(Number.isFinite(x)&&Number.isFinite(y));assert.ok(x>=90&&x<=1170);
  }
 }
 assert.equal(run('Object.keys(state.models).length'),0);assert.equal(run('Object.keys(state.regions).length'),0);
 assert.equal(JSON.stringify(run('DATA')),before);assert.deepEqual(dense.errors,[]);dense.dom.window.close();
 // A correspondence cannot be inferred through a comparison row with no hit.
 const gap=fixture(1);gap.regions[2].genome_id='Species_c_GCF_3.1';
 gap.models=gap.models.map(m=>m.region_id==='r2'?{...m,genome_id:gap.regions[2].genome_id,gene_group:m.gene_group.replace(/^a_/,'c_')}:m);
 gap.neighbor_pairs=[{...gap.neighbor_pairs.find(p=>p.gene_a==='a_1'),gene_b:'c_1',model_b:'r2_1_annotation'}];
 gap.source_neighbor_pairs.tiberius=gap.neighbor_pairs.map(p=>({...p,model_a:p.model_a.replace(/annotation$/,'tiberius'),model_b:p.model_b.replace(/annotation$/,'tiberius')}));
 const g=open(gap);g.change('other','all');assert.equal(g.d.querySelectorAll('#links .neighbor-link').length,0);assert.deepEqual(g.errors,[]);g.dom.window.close();
 // A failed/empty region still offers a review handoff; it does not invent an
 // overlapping model or silently label the gene absent.
 const empty=fixture(0);empty.regions=empty.regions.slice(0,1);empty.models=[];empty.neighbor_pairs=[];empty.errors=[{stage:'context_prediction',error:'deliberate fixture failure'}];
 const e=open(empty);assert.equal(e.d.querySelectorAll('[data-match-type="empty-anchor"]').length,1);assert.match(e.d.getElementById('models').textContent,/does not establish absence/);
 e.run('window.__downloads=[];download=(name,text)=>window.__downloads.push({name,text})');
 e.change('regionDecision','unresolved');e.d.getElementById('regionNotes').value='Isolated failure control';e.d.getElementById('regionNotes').dispatchEvent(new e.w.Event('input'));e.d.getElementById('saveRegion').click();e.d.getElementById('export').click();
 assert.equal(e.run('window.__downloads.length'),1);assert.deepEqual(e.errors,[]);e.dom.window.close();
 const context=fixture(0);context.annotation_context_available=true;
 context.annotation_features=[{feature_id:'rna',region_id:'r0',genome_id:context.regions[0].genome_id,original_gene:'rna',name:'TEST_LNCRNA',kind:'gene',biotype:'lncRNA',contig:'chr',start:499000,end:499300,strand:'+'}];
 const c=open(context);c.change('scale','models');assert.equal(c.d.querySelectorAll('#links .gene-arrow[data-row="0"]').length,2);
 assert.equal(c.d.querySelector('#links .gene-count[data-row="0"]').dataset.annotationGenes,'2');
 assert.equal(c.d.querySelectorAll('#links .source-arrow[data-context-only="true"]').length,1);
 assert.equal(c.d.querySelectorAll('#models .context-feature select').length,0,'A noncoding GFF feature must never become a selectable protein model');
 c.change('modelSource','tiberius');assert.equal(c.d.querySelectorAll('#links .gene-arrow[data-row="0"]').length,1);
 assert.equal(c.run('Object.keys(state.models).length'),0);assert.deepEqual(c.errors,[]);c.dom.window.close();
 // Exact duplicates collapse visually; overlap, equal CDS strings, phase changes,
 // different strands and differences outside the clipped window stay distinct.
 const compact=fixture(0);compact.regions=compact.regions.slice(0,1);compact.models=compact.models.filter(m=>m.region_id==='r0');compact.neighbor_pairs=[];
 const q=open(compact);assert.equal(q.d.querySelector('.cluster-gene-label').dataset.sourceBadge,'A+T');assert.equal(q.d.querySelectorAll('.cluster-model-span').length,1);q.change('scale','models');
 assert.equal(q.d.querySelectorAll('#links .ruler-gene-label').length,1);
 assert.equal(q.d.querySelectorAll('#links .transcript-model').length,1);
 assert.equal(q.d.querySelector('.model-source-badge').childNodes[0].textContent,'A+T');
 assert.equal(q.d.querySelector('.transcript-model').dataset.modelIds.split(' ').length,2);
 q.change('modelSource','annotation');assert.equal(q.d.querySelector('.model-source-badge').childNodes[0].textContent,'A');
 q.change('modelSource','tiberius');assert.equal(q.d.querySelector('.model-source-badge').childNodes[0].textContent,'T');
 q.dom.window.close();
 for(const kind of ['coordinates','phase','strand','cds_hash','outside_window','different_gene']){
  const data=JSON.parse(JSON.stringify(compact)),a=data.models[0],t=data.models[1];
  if(kind==='coordinates'){t.start+=10;t.cds_blocks[0][0]+=10}
  if(kind==='phase')t.cds_blocks[0][2]=1;
  if(kind==='strand')t.strand=t.strand==='+'?'-':'+';
  if(kind==='cds_hash'){a.cds_sha256='aaa';t.cds_sha256='bbb'}
  if(kind==='outside_window'){
   a.end=1000100;t.end=1000200;a.cds_blocks.push([999900,1000100,0]);t.cds_blocks.push([999900,1000200,0]);
  }
  if(kind==='different_gene')t.gene_group='different_physical_gene';
  const v=open(data);if(kind!=='different_gene')assert.match(v.d.querySelector('.gene-arrow title').textContent,/A\/T/);assert.equal(v.d.querySelectorAll('.cluster-model-span').length,2);v.change('scale','models');assert.equal(v.d.querySelectorAll('#links .transcript-model').length,2,kind+' must not collapse');assert.deepEqual(v.errors,[]);v.dom.window.close();
 }
 // Source-only target models must stay visible and filter emptiness must be explicit.
 const onlyT=fixture(0);onlyT.regions=onlyT.regions.slice(0,1);onlyT.models=onlyT.models.filter(m=>m.region_id==='r0'&&m.source==='tiberius');onlyT.neighbor_pairs=[];
 const f=open(onlyT,'?locus=r0');
 assert.equal(f.d.getElementById('region').value,'r0');assert.equal(f.d.querySelectorAll('.cluster-model-span').length,1);
 const targetShape=f.d.querySelector('.cluster-model-span');assert.ok(Number(targetShape.getAttribute('stroke-width'))>=1.8);assert.ok(Number(targetShape.getAttribute('fill-opacity'))>0);
 f.change('modelSource','annotation');assert.equal(f.d.querySelectorAll('.cluster-model-span').length,0);assert.equal(f.d.getElementById('diagramStatus').hidden,false);assert.match(f.d.getElementById('diagramStatus').textContent,/filter hides available Tiberius/);assert.match(f.d.querySelector('[data-match-type="empty-anchor"]').textContent,/filter hides available/);
 f.d.getElementById('diagramStatus').querySelector('button').click();assert.equal(f.d.getElementById('modelSource').value,'all');assert.equal(f.d.querySelectorAll('.cluster-model-span').length,1);assert.equal(f.d.getElementById('diagramStatus').hidden,true);assert.equal(f.run('Object.keys(state.models).length'),0);f.dom.window.close();
 const linked=open(fixture(0),'?locus=r1');assert.equal(linked.d.getElementById('region').value,'r1');assert.equal(linked.run('Object.keys(state.models).length'),0);linked.dom.window.close();
 const invalid=open(fixture(0),'?locus=missing');assert.equal(invalid.d.getElementById('region').value,'r0');assert.deepEqual(invalid.errors,[]);invalid.dom.window.close();
 const result={status:'passed',source_filter_empty_state_explicit:true,tiberius_only_contrast:true,locus_deep_link_without_decisions:true,exact_structure_collapse_only:true,source_badges_and_filters:true,synthetic_gene_groups_per_region:45,default_displayed_groups_per_region:21,model_alternatives:270,dense_both_scales:true,hard_neighbor_display_limit:true,no_bridging_missing_correspondence:true,empty_failed_region_handoff:true,no_production_decisions_modified:true,elapsed_ms:performance.now()-start};
 if(process.argv[2])fs.writeFileSync(process.argv[2],JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));
})().catch(e=>{console.error(e);process.exitCode=1});
