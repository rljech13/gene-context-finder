'use strict';
const fs=require('node:fs'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const {JSDOM,VirtualConsole}=require('jsdom');const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
const html=fs.readFileSync(process.argv[2],'utf8');
// Parse the exact embedded JSON directly: jsdom's HTML tokenizer otherwise
// creates millions of character tokens for the largest offline data literal.
const start=html.indexOf('const DATA=hydrateReport(')+'const DATA=hydrateReport('.length,end=html.indexOf(');\nconst $=',start);assert.ok(start>25&&end>start);
const exactPayload=JSON.parse(html.slice(start,end)),runtimeHTML=html.slice(0,start)+'__savedBlastPayload'+html.slice(end);
const dom=new JSDOM(runtimeHTML,{runScripts:'dangerously',url:'http://localhost/',virtualConsole:vc,beforeParse(w){w.__savedBlastPayload=exactPayload}}),w=dom.window,d=w.document,run=s=>w.eval(s);
const dataFingerprint=()=>crypto.createHash('sha256').update(run(`(()=>{const seen=new WeakMap();let next=0;function visit(v){if(!v||typeof v!=='object')return v;if(seen.has(v))return {$runtime_ref:seen.get(v)};seen.set(v,next++);return Array.isArray(v)?v.map(visit):Object.fromEntries(Object.entries(v).map(([k,x])=>[k,visit(x)]))}return JSON.stringify(visit(DATA))})()`)).digest('hex');
assert.deepEqual(errors,[]);const before=dataFingerprint(),draft=run('JSON.stringify(draft())');
const counts=run(`(()=>{let hsps=0;for(const l of DATA.loci)for(const c of l.chains)for(const h of c.hsps){const a=h.alignment||DATA.blast_alignments?.[l.genome_id]?.[h.hsp_id];if(!a)throw Error('Missing HSP trace');if(a.qseq.length!==h.alignment_length||a.sseq.length!==h.alignment_length)throw Error('Bad columns');hsps++}return {loci:DATA.loci.length,hsps}})()`);
if(!d.querySelector('.candidate'))run("$('genome').value=DATA.genomes.find(g=>groups.get(g.genome_id).length).genome_id;render()");
const toggle=[...d.querySelectorAll('summary')].find(s=>s.textContent==='BLAST alignment — query vs candidate');assert.ok(toggle);toggle.parentElement.open=true;toggle.parentElement.dispatchEvent(new w.Event('toggle'));
const box=d.querySelector('.blast-alignment');assert.ok(box);assert.ok(box.querySelector('.blast-lines'));assert.ok(box.querySelector('.blast-match'));
const query=box.querySelector('[aria-label="BLAST query"]'),product=box.querySelector('[aria-label="BLAST candidate product"]');
let switches=0;for(const o of [...query.options]){query.value=o.value;query.dispatchEvent(new w.Event('change'));assert.ok(box.querySelector('.blast-lines'));if(product.options.length>1){product.selectedIndex=1;product.dispatchEvent(new w.Event('change'));assert.ok(box.querySelector('.blast-lines'));}switches++}
// A real gapped trace must retain dashes, exact query residues and endpoint positions.
const gap=run(`(()=>{for(const l of DATA.loci)for(const c of l.chains)for(const h of c.hsps){const a=h.alignment||DATA.blast_alignments?.[l.genome_id]?.[h.hsp_id];if(a.qseq.includes('-')||a.sseq.includes('-'))return {l,c,h,a}}return null})()`);
if(gap){w.__l=gap.l;const b=run('blastAlignment(__l)');assert.ok(b);}
assert.equal(dataFingerprint(),before);assert.equal(run('JSON.stringify(draft())'),draft);assert.deepEqual(errors,[]);
const result={status:'passed',sha256:crypto.createHash('sha256').update(html).digest('hex'),...counts,query_switches:switches,product_selector:true,gapped_trace_present:!!gap,scientific_data_unchanged:true,decisions_unchanged:true,exact_embedded_payload_injected:true,errors};if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2));console.log(JSON.stringify(result));dom.window.close();
