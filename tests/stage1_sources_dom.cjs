/* One source-context presentation contract for any gene panel. */
const fs=require('node:fs'),assert=require('node:assert/strict'),{JSDOM,VirtualConsole}=require('jsdom');
const results=[];
for(const file of process.argv.slice(2)){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(file,'utf8'),{runScripts:'dangerously',url:'http://sources.invalid/',virtualConsole:vc,beforeParse(w){w.HTMLElement.prototype.scrollIntoView=function(){}}}),w=dom.window,d=w.document,run=s=>w.eval(s),data=run('DATA'),before=run('JSON.stringify(DATA)');
 assert.deepEqual(errors,[]);assert.ok(data.queries.length);const origins=data.early_synteny.source_resolution,missing=origins.filter(r=>!r.sources.length),bound=origins.filter(r=>r.context_binding_only);
 assert.equal(d.getElementById('sourceWarnings').hidden,!(missing.length||bound.length));
 let views=0,controls=0;
 for(const genome of data.genomes){
  d.getElementById('genome').value=genome.genome_id;d.getElementById('genome').dispatchEvent(new w.Event('change'));
  assert.ok(d.querySelectorAll('.candidate').length<=10);const l=data.loci.find(l=>l.genome_id===genome.genome_id&&l.early_synteny.comparisons.length);if(!l){const unsupported=data.loci.find(l=>l.genome_id===genome.genome_id);if(unsupported){w.__l=unsupported;const text=run('synteny(__l).textContent');if(unsupported.early_synteny.status==='annotation_unavailable')assert.ok(text.includes('candidate assembly has no supplied annotation'))}continue;}
  w.__l=l;const view=run('synteny(__l)'),toggle=view.querySelector('input[aria-label="Show self comparison"]');assert.ok(toggle);assert.equal(toggle.checked,false);
  const headings=()=>[...view.querySelectorAll('svg text')].map(n=>n.textContent).filter(t=>/^(Candidate|Reference|Same assembly|Self-control):/.test(t));
  const normal=headings();assert.equal(normal.filter(t=>t.startsWith('Self-control:')).length,0);
  const self=l.early_synteny.comparisons.filter(c=>{w.__r=run(`refsById.get(${JSON.stringify(c.reference_id)})`);return run('comparisonRole(__l,__r)')==='self-control'}).length;
  const cross=l.early_synteny.comparisons.filter(c=>{w.__r=run(`refsById.get(${JSON.stringify(c.reference_id)})`);return run('comparisonRole(__l,__r)')==='cross-species'}).length;
  assert.equal(normal.filter(t=>t.startsWith('Reference:')).length,cross);
  toggle.checked=true;toggle.dispatchEvent(new w.Event('change'));assert.equal(headings().filter(t=>t.startsWith('Self-control:')).length,self);controls+=self;views++;
 }
 assert.equal(run('JSON.stringify(DATA)'),before);assert.deepEqual(errors,[]);
 results.push({file,queries:data.queries.length,source_loci:data.early_synteny.references.length,resolved_variants:origins.length-missing.length,missing_variants:missing.length,source_sequence_disagreements:bound.length,tested_genomes:views,controls_hidden_by_default:controls,scientific_payload_unchanged:true,errors});dom.window.close();
}
console.log(JSON.stringify(results,null,2));
