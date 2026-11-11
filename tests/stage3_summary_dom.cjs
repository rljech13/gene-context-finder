/* Compact report: all-by-all alignments followed by the exact selected set. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{JSDOM,VirtualConsole}=require('jsdom');
const root=path.resolve(process.argv[2]),errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(path.join(root,'report.html'),'utf8'),{runScripts:'dangerously',url:'http://isolated-sheet.invalid/report.html',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){};w.TextEncoder=global.TextEncoder}});
const w=dom.window,d=w.document,run=s=>w.eval(s),before=run('JSON.stringify(DATA)'),change=(id,v)=>{const n=d.getElementById(id);n.value=v;n.dispatchEvent(new w.Event('change'))};
const models=run('selected()'),proteinRefs=run("refs('protein')");
assert.deepEqual(errors,[]);assert.ok(models.length);
assert.equal(d.getElementById('referenceTable'),null,'The initial comparison table was removed');
assert.equal(d.getElementById('comparisonTitle'),null);
assert.equal(d.getElementById('modelSearch'),null);
assert.equal(d.querySelectorAll('details[open]').length,1);assert.equal(d.getElementById('alignmentPanel').open,true);
assert.equal(d.getElementById('groupLabel').hidden,run('DATA.groups.length')<2);
const aligned=run('DATA.alignments[group()].protein'),ids=run("alignmentIds('protein')"),length=Object.values(aligned)[0].length;
assert.equal(d.getElementById('hideEmptyColumns'),null);
const kept=Array.from({length},(_,i)=>i).filter(i=>ids.some(id=>aligned[id][i]!=='-')).length;
assert.equal(d.querySelectorAll('#alignment .aln-row').length,ids.length);
for(const row of d.querySelectorAll('#alignment .aln-row'))assert.equal(row.querySelectorAll('.res').length,kept);
assert.equal(d.querySelectorAll('#allPairMatrix tbody tr').length,ids.length);
assert.equal(d.querySelectorAll('#allPairMatrix tbody td').length,ids.length**2);
for(const [a,b] of [[models[0].model_id,models[2].model_id],[models[2].model_id,models[0].model_id],[models[0].model_id,proteinRefs[0].id],[proteinRefs[0].id,proteinRefs[1].id]]){
 run(`openPair(${JSON.stringify(a)},${JSON.stringify(b)},'protein')`);
 assert.equal(d.getElementById('alignmentModel').value,a);assert.equal(d.getElementById('alignmentReference').value,b);
 assert.equal(d.querySelectorAll('#alignment .aln-row').length,2);assert.ok(d.querySelector('#pairDownload a'));
 const records=[...d.querySelectorAll('#alignment .aln-row code')].map(n=>n.textContent);
 assert.equal(records[0].replaceAll('-',''),aligned[a].replaceAll('-',''));assert.equal(records[1].replaceAll('-',''),aligned[b].replaceAll('-',''));
 for(let i=0;i<records[0].length;i++)assert.notEqual(records[0][i]+records[1][i],'--');
}
const matrixSelf=d.querySelector('#allPairMatrix tbody tr td button');assert.equal(matrixSelf.textContent,'100.0%');matrixSelf.click();assert.equal(d.querySelectorAll('#alignment .aln-row').length,1);
change('alignmentKind','cds');change('alignmentModel','all');const cds=run('DATA.alignments[group()].cds'),cdsIds=run("alignmentIds('cds')"),cdsLength=Object.values(cds)[0].length;
const cdsKept=Array.from({length:cdsLength},(_,i)=>i).filter(i=>cdsIds.some(id=>cds[id][i]!=='-')).length;
for(const row of d.querySelectorAll('#alignment .aln-row'))assert.equal(row.querySelectorAll('.res').length,cdsKept);
assert.equal(d.querySelectorAll('#allPairMatrix tbody td').length,cdsIds.length**2);
change('allMatrixMetric','coverage');assert.ok(!/NaN|undefined|Infinity/.test(d.getElementById('allPairMatrix').textContent));
const final=run('finalRows()');assert.equal(final.primary.length+final.cases.length,models.length);
assert.equal(d.querySelectorAll('#decisionTable tbody tr').length,final.primary.length);
assert.equal(d.querySelectorAll('#specialCaseTable tbody tr').length,final.cases.length);
assert.equal(new Set(final.primary.map(r=>r.group+'|'+r.genome_id)).size,final.primary.length);
assert.deepEqual([...d.querySelectorAll('#finalChoices [data-retained-model-id]')].map(n=>n.dataset.retainedModelId).sort(),Array.from(models,m=>m.model_id).sort());
assert.ok(d.getElementById('alignmentPanel').compareDocumentPosition(d.getElementById('finalChoices'))&w.Node.DOCUMENT_POSITION_FOLLOWING);
for(const r of [...final.primary,...final.cases]){assert.ok(r.selection_reason);assert.equal(r.decision_origin,'suggestion')}
const ds=d.getElementById('decisionSearch');ds.value=models[0].species;ds.dispatchEvent(new w.Event('input'));
assert.equal(d.querySelectorAll('#finalChoices tbody tr').length,models.filter(m=>m.species===models[0].species).length);
ds.value='no-matching-organism';ds.dispatchEvent(new w.Event('input'));assert.equal(d.querySelectorAll('#finalChoices tbody tr').length,0);ds.value='';ds.dispatchEvent(new w.Event('input'));
const bucket=run('finalSpec.buckets.find(b=>b.model_ids.length>1)');assert.ok(bucket);const original=run(`finalState.choices[${JSON.stringify(bucket.key)}].primary`),replacement=bucket.model_ids.find(id=>id!==original);
run(`editPrimary(${JSON.stringify(replacement)})`);assert.equal(d.getElementById('finalModel').value,replacement);assert.equal(d.getElementById('finalReason').value,'');d.getElementById('finalSave').click();assert.equal(run(`finalState.choices[${JSON.stringify(bucket.key)}].primary`),original);d.getElementById('finalReason').value='Reviewed alternative after locus inspection';d.getElementById('finalSave').click();
assert.equal(run(`finalState.choices[${JSON.stringify(bucket.key)}].primary`),replacement);assert.ok(d.querySelector(`#decisionTable tr[data-retained-model-id="${replacement}"]`));assert.ok(d.querySelector(`#specialCaseTable tr[data-retained-model-id="${original}"]`));
assert.equal(run(`JSON.parse(localStorage.getItem(finalKey)).choices[${JSON.stringify(bucket.key)}].reviewed`),true);
run(`editPrimary(${JSON.stringify(replacement)})`);change('finalModel','');assert.equal(d.getElementById('finalReason').value,'');d.getElementById('finalReason').value='Copy identity unresolved';d.getElementById('finalSave').click();assert.equal(run(`finalState.choices[${JSON.stringify(bucket.key)}].primary`),'');
for(const id of bucket.model_ids)assert.ok(d.querySelector(`#specialCaseTable tr[data-retained-model-id="${id}"]`));
run(`editCase(${JSON.stringify(original)})`);d.getElementById('finalReason').value='A retained separate copy';d.getElementById('finalSave').click();assert.equal(run(`finalRows().cases.find(r=>r.model_id===${JSON.stringify(original)}).selection_reason`),'A retained separate copy');
for(const mutate of [s=>s.set_id='wrong',s=>s.review_id='wrong',s=>s.choices[bucket.key].primary='reference_q0001',s=>s.choices[bucket.key].reason='',s=>s.choices[bucket.key].reviewed='true']){const invalid=JSON.parse(run('JSON.stringify(finalState)'));mutate(invalid);assert.throws(()=>run(`validateFinal(${JSON.stringify(invalid)})`))}
const files=run('finalFiles()');for(const key of ['primary_targets.tsv','special_cases.tsv','stage3.final.decisions.json','primary_targets.cds.fna','primary_targets.proteins.faa','special_cases.cds.fna','special_cases.proteins.faa'])assert.ok(key in files);
for(const[name,text]of Object.entries(files))if(name.endsWith('.fna')||name.endsWith('.faa')){const kind=name.endsWith('.fna')?'cds':'protein',rows=name.startsWith('primary')?run('finalRows().primary'):run('finalRows().cases');const records=new Map(text.split('>').slice(1).map(r=>{const[h,...seq]=r.trim().split('\n');return[h,seq.join('')]}));assert.equal(records.size,rows.length);for(const r of rows)assert.equal(records.get(r.model_id),run(`DATA.finalization.raw_sequences.${kind}[${JSON.stringify(r.model_id)}]`))}
assert.ok(run('zipFiles(finalFiles()).size')>1000);
if(process.argv[3]){fs.writeFileSync(process.argv[3]+'.edited.json',run('JSON.stringify(finalState,null,2)'));fs.writeFileSync(process.argv[3]+'.browser-files.json',JSON.stringify(files))}
let links=0;for(const n of d.querySelectorAll('a[href]')){const url=n.getAttribute('href');if(!url||url.startsWith('#')||/^https?:/.test(url))continue;assert.ok(fs.existsSync(path.join(root,decodeURIComponent(url.split(/[?#]/)[0]))),url);links++}
assert.equal(d.querySelectorAll('script[src],link[rel=stylesheet],iframe').length,0);
assert.equal(run('JSON.stringify(DATA)'),before);assert.deepEqual(errors,[]);
const result={status:'passed',selected_models:models.length,protein_references:proteinRefs.length,initial_reference_table_removed:true,all_by_all_pairs:true,pair_alignment_downloads_verified:true,self_comparisons_visible:true,empty_columns_always_hidden:true,final_choices_after_alignments:true,primary_targets:final.primary.length,separate_cases:final.cases.length,one_primary_per_organism:true,editable_assignment:true,exact_partition_exports:true,portable_decisions:true,sequence_search:true,offline_local_links_checked:links,scientific_payload_unchanged:true,errors};
if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));dom.window.close();
