/* Offline interaction regression: NODE_PATH=/tmp/gcf-ui/node_modules node tests/stage1_groups_dom.cjs REPORT [OUTPUT] */
'use strict';
const fs=require('node:fs'),assert=require('node:assert/strict');
const {JSDOM,VirtualConsole}=require('jsdom');
const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const dom=new JSDOM(fs.readFileSync(process.argv[2],'utf8'),{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc});
const w=dom.window,d=w.document,run=s=>w.eval(s),plain=x=>JSON.parse(JSON.stringify(x));
assert.deepEqual(errors,[]);
const before=run('JSON.stringify(DATA)'),legacy=plain(run('draft()'));
const gs=plain(run('DATA.genomes.filter(g=>groups.get(g.genome_id).length>=2).slice(0,2).map(g=>({id:g.genome_id,loci:groups.get(g.genome_id).slice(0,2).map(l=>l.locus_id)}))'));
assert.equal(gs.length,2);
for(const p of legacy.pairs){p.selected=gs.find(g=>g.id===p.genome_id)?.loci||[];p.action=p.selected.length?'inspect':'unresolved';p.reason=p.selected.length?'compare_copies':'insufficient_evidence';p.committed=true;}
w.__draft=legacy;run('loadDraft(__draft);render()');
// Create named sets through the real form.
for(const name of ['Family A','Family B','Family C']){d.getElementById('analysisGroupName').value=name;d.getElementById('addAnalysisGroup').click();}
assert.equal(run('analysisGroups.length'),3);
const assignments=[[gs[0].loci[0],'group-1'],[gs[1].loci[0],'group-1'],[gs[0].loci[1],'group-2'],[gs[1].loci[1],'group-3']];
for(const [id,gid] of assignments){const select=d.querySelector(`select[aria-label="Analysis group for ${id}"]`);assert.ok(select);select.value=gid;select.dispatchEvent(new w.Event('change'));}
const full=plain(run('draft()'));assert.equal(full.schema,'stage1_pair_draft_v2');assert.equal(full.analysis_groups.length,3);
for(const p of legacy.pairs){const updated=full.pairs.find(x=>x.genome_id===p.genome_id);assert.deepEqual(updated.selected,p.selected);assert.equal(updated.committed,true);}
assert.deepEqual(plain(run('analysisGroupErrors()')),[]);
const groupExports=[];
for(const gid of ['group-1','group-2','group-3']){w.__gid=gid;const value=plain(run('groupDraft(__gid)'));assert.equal(value.active_analysis_group,gid);assert.equal(value.pairs.length,legacy.pairs.length);groupExports.push(value);}
// Download button uses the scoped document, not the global draft.
let clicked=null;w.download=(name,text)=>{clicked={name,value:JSON.parse(text)}};
d.querySelector('button[aria-label="Download group Family A"]').click();assert.equal(clicked.value.active_analysis_group,'group-1');
w.__draft=clicked.value;run('loadDraft(__draft)');d.getElementById('draft').click();assert.equal(clicked.value.active_analysis_group,undefined);assert.equal(clicked.value.analysis_groups.length,3);
// Reload + TSV retain membership and all old selected target IDs.
w.__draft=full;run('loadDraft(__draft)');
const rows=plain(run('exportRows()'));assert.equal(rows.filter(r=>r.analysis_group_id).length,4);
w.__rows=rows;const tsv=run("[columns.join('\\t'),...__rows.map(r=>columns.map(k=>String(r[k]??'').replace(/[\\t\\r\\n]/g,' ')).join('\\t'))].join('\\n')+'\\n'");
w.__tsv=tsv;run('loadTSV(__tsv)');assert.deepEqual(plain(run('draft().analysis_groups')),full.analysis_groups);
// Failed imports are atomic; assignment never approves an unsaved pair.
const bad=structuredClone(full);bad.analysis_groups[1].members.push(bad.analysis_groups[0].members[0]);w.__draft=bad;
const saved=run('JSON.stringify(draft())');assert.throws(()=>run('loadDraft(__draft)'),/multiply assigned/);assert.equal(run('JSON.stringify(draft())'),saved);
w.__gid=gs[0].id;run('state.get(__gid).committed=false');assert.throws(()=>run("groupDraft('group-1')"),/not been saved/);
w.__lid=gs[0].loci[0];run("assignAnalysisGroup(__lid,'group-2')");assert.equal(run('state.get(__gid).committed'),false);
// Legacy draft resets group metadata without changing the old choices.
w.__draft=legacy;run('loadDraft(__draft)');assert.equal(run('analysisGroups.length'),0);assert.equal(run('draft().schema'),'stage1_pair_draft_v1');
assert.equal(run('JSON.stringify(DATA)'),before);assert.deepEqual(errors,[]);
const receipt={status:'passed',groups:3,selected_loci:4,same_assembly_multiple_groups:true,download_scopes_group:true,tsv_round_trip:true,legacy_compatible:true,uncommitted_rejected:true,invalid_import_atomic:true,scientific_data_unchanged:true,errors};
if(process.argv[3]){fs.writeFileSync(process.argv[3],JSON.stringify(receipt,null,2)+'\n');fs.writeFileSync(process.argv[3].replace(/\.json$/,'.drafts.json'),JSON.stringify(groupExports,null,2)+'\n');}
console.log(JSON.stringify(receipt));dom.window.close();
