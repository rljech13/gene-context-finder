"""Offline locus review: filter the display, inspect tracks, export decisions."""
import json


def write_report(root, data):
    payload = json.dumps(data).replace('<', '\\u003c').replace('&', '\\u0026')
    (root/'report.html').write_text(TEMPLATE.replace('__DATA__', payload))


TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Stage 1 — candidate loci</title>
<style>
:root{color-scheme:light;--ink:#173044;--muted:#536575;--line:#cbd5df;--paper:#f5f7fa;--blue:#245cc4}
*{box-sizing:border-box}body{font:15px/1.5 system-ui,sans-serif;color:var(--ink);background:var(--paper);margin:0}
main{max-width:1450px;margin:auto;padding:28px}h1{margin:0 0 8px;font-size:30px}h2{font-size:20px;margin:0 0 12px}
p{margin:8px 0 14px}.muted{color:var(--muted)}.notice{background:#fff4d7;padding:14px;border-radius:8px;margin:18px 0}
.toolbar,.card,.summary{background:white;border:1px solid var(--line);border-radius:10px;padding:18px;margin:16px 0}
.toolbar{position:sticky;top:0;z-index:2;box-shadow:0 4px 12px #17304412}.controls{display:flex;gap:15px;flex-wrap:wrap;align-items:end}
label{display:flex;flex-direction:column;gap:4px;font-size:13px}input,select,button{font:inherit;padding:7px 9px;border:1px solid #a8bacb;border-radius:5px}
input[type=number]{width:120px}input[type=search]{width:240px}button{cursor:pointer;background:var(--blue);color:white;border:0}
button.secondary{background:#e8eef9;color:#183c77}label.check{flex-direction:row;align-items:center;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--muted)}code{overflow-wrap:anywhere}.scroll{overflow-x:auto}.tracks{width:100%;min-width:850px;height:auto}
.badge{display:inline-block;padding:2px 8px;border-radius:12px;background:#e7effe;margin:2px;font-size:12px}
.bad{background:#fff0cc}.good{background:#dff3e8}.decision{display:flex;flex-wrap:wrap;gap:12px;margin-top:14px;align-items:end}
.decision .notes{flex:1;min-width:220px}.decision .notes input{width:100%}.anchor input{width:230px}.hidden{display:none}
.error{color:#a02525}a{color:#245cc4}details{margin-top:12px}summary{cursor:pointer;font-weight:600}.legend{font-size:12px;color:var(--muted)}
@media(max-width:700px){main{padding:12px}.toolbar{position:static}.controls{gap:8px}}
</style></head><body><main>
<h1 id="title">Stage 1 — candidate loci</h1>
<p>Review assembled loci, their reference support and proposed intervals. No locus is approved automatically.</p>
<div class="notice">Coordinates in this report and exported TSV/BED are <strong>zero-based, half-open</strong>.
BLAST chains are locus proposals; local miniprot models are predictions. Similarity and filter status do not establish orthology.
The tracks are alignments of the supplied references, not a general gene annotation.</div>
<div class="summary"><h2>Genome coverage</h2><div id="genomes" class="scroll"></div>
<p class="muted">No candidate means no retained candidate under these search/assembly settings; it does not establish gene absence.</p></div>
<div class="toolbar"><div class="controls">
<label>Find genome / locus / contig<input id="find" type="search" placeholder="All genomes"></label>
<label>Query coverage ≥<input id="coverage" type="number" min="0" max="1" step="0.05"></label>
<label>Identity ≥<input id="identity" type="number" min="0" max="1" step="0.05"></label>
<label>E-value ≤<input id="evalue" type="number" min="0" step="any"></label>
<label class="check"><input id="passing" type="checkbox">Show passing loci only</label>
<button id="export">Download decisions TSV</button>
<label>Reload saved decisions<input id="import" type="file" accept=".tsv,.txt"></label>
</div><p id="counts" class="muted"></p>
<p class="legend">Filters change this display only. Export includes decisions for <strong>all</strong> loci and genomes without candidates,
including hidden rows. Use <code>stage1-filter</code> to save a new filtered result. Save/download decisions before closing.</p></div>
<div id="loci"></div>
<details class="summary"><summary>Saved settings and files</summary><pre id="settings"></pre>
<p><a href="loci.tsv">All loci TSV</a> · <a href="passing_loci.tsv">Saved passing loci TSV</a> ·
<a href="query_support.tsv">Per-reference support</a> · <a href="loci.fna">Candidate DNA FASTA</a> ·
<a href="loci.decisions.tsv">Blank decision template</a> · <a href="manifest.json">Provenance and commands</a></p></details>
</main><script>
const DATA=__DATA__;
const $=id=>document.getElementById(id), queries=Object.fromEntries(DATA.queries.map(q=>[q.id,q]));
const columns=['review_id','genome_id','locus_id','decision','anchor_intervals','window_start','window_end','notes'];
const entries=[], percent=x=>(100*x).toFixed(1)+'%';
function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n}
function table(headers,rows){const t=el('table'),head=el('tr');headers.forEach(h=>head.append(el('th',h)));t.append(head);
 for(const row of rows){const tr=el('tr');row.forEach(x=>tr.append(el('td',String(x))));t.append(tr)}return t}
function svgEl(tag,attrs,text){const n=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n}
function tracks(l,full=false){
 let chains=l.chains,models=l.models||[];
 if(!full){const best=new Map();for(const c of chains){const member=queries[c.query].members[0],key=(member.species||c.query)+' '+(member.gene_symbol||'');
 if(!best.has(key)||c.coverage>best.get(key).coverage)best.set(key,c)}chains=[...best.values()];models=models.slice(0,1)}
 const items=[...chains.map(c=>({kind:'BLAST',query:c.query,blocks:c.blocks,coverage:c.coverage,identity:c.identity})),
 ...models.map(m=>({kind:'CDS?',query:m.query,blocks:m.cds.map(c=>c.slice(0,2)),coverage:m.coverage,identity:m.identity}))];
 const lower=Math.min(l.start,...models.map(m=>m.start)),upper=Math.max(l.end,...models.map(m=>m.end));
 const pad=Math.max(30,Math.round((upper-lower)*.04)),lo=Math.max(0,lower-pad),hi=Math.min(l.contig_length,upper+pad);
 const left=300,right=1160,scale=x=>left+(x-lo)/(hi-lo)*(right-left),svg=svgEl('svg',{viewBox:`0 0 1200 ${65+items.length*32}`,class:'tracks',role:'img','aria-label':'Reference alignment and predicted CDS tracks'});
 svg.append(svgEl('line',{x1:left,y1:22,x2:right,y2:22,stroke:'#8399ad'}));
 svg.append(svgEl('text',{x:left,y:14,'font-size':12},lo.toLocaleString()));svg.append(svgEl('text',{x:right,y:14,'text-anchor':'end','font-size':12},hi.toLocaleString()));
 items.forEach((item,i)=>{const y=50+i*32,q=queries[item.query],members=q.members;
 const label=`${item.kind} · ${members[0].species||item.query} · ${members[0].gene_symbol||item.query}`;
 const text=svgEl('text',{x:0,y:y+4,'font-size':12},label);text.append(svgEl('title',{},members.map(m=>m.query).join('\n')));svg.append(text);
 const a=Math.min(...item.blocks.map(b=>b[0])),b=Math.max(...item.blocks.map(b=>b[1]));
 svg.append(svgEl('line',{x1:scale(a),y1:y,x2:scale(b),y2:y,stroke:'#9eacb9'}));
 item.blocks.forEach(block=>{const rect=svgEl('rect',{x:scale(block[0]),y:y-6,width:Math.max(1,scale(block[1])-scale(block[0])),height:12,fill:item.kind==='BLAST'?'#3266c8':'#3a9377'});
 rect.append(svgEl('title',{},`${block[0]}:${block[1]} | query coverage ${percent(item.coverage)}, identity ${percent(item.identity)}`));svg.append(rect)});
 });return svg;
}
function persist(){try{localStorage.setItem('gcf-stage1-'+DATA.review_id,JSON.stringify(entries.map(e=>e.row)))}catch(_){}}
function input(label,field,row,cls){const wrap=el('label',label,cls),n=el('input');n.value=row[field]||'';n.addEventListener('input',()=>{row[field]=n.value;persist()});wrap.append(n);return {wrap,n}}
function decision(card,row,l){const box=el('div',undefined,'decision'),choice=el('label','Decision'),select=el('select');
 [['','Pending'],['inspect','Inspect further'],['exclude','Exclude'],['unresolved','Unresolved']].forEach(([value,label])=>{const o=el('option',label);o.value=value;select.append(o)});
 select.addEventListener('change',()=>{row.decision=select.value;persist()});choice.append(select);box.append(choice);
 let anchor;
 if(l){anchor=input('Anchor start:end (semicolon to split)', 'anchor_intervals',row,'anchor');box.append(anchor.wrap);
 box.append(input('Custom window start','window_start',row).wrap,input('Custom window end','window_end',row).wrap)}
 const notes=input('Notes (required for exclusions, unresolved or overrides)','notes',row,'notes');box.append(notes.wrap);card.append(box);
 const entry={row,card,select,anchor:anchor?.n,notes:notes.n,locus:l};entries.push(entry);return entry;
}
$('title').textContent=`${DATA.config.gene} — stage 1 candidate loci`;
$('coverage').value=DATA.filters.min_query_coverage;$('identity').value=DATA.filters.min_identity;$('evalue').value=DATA.filters.max_evalue;
$('settings').textContent=JSON.stringify({search_and_assembly:DATA.config,saved_review_filters:DATA.filters,review_id:DATA.review_id},null,2);
$('genomes').append(table(['Genome','Search status','Raw HSPs','Candidate loci','Expected control overlap / error'],DATA.genomes.map(g=>[
 g.genome_id,g.status,g.raw_hsps||0,DATA.loci.filter(l=>l.genome_id===g.genome_id).length,g.error||(g.refinement_errors?.length?g.refinement_errors.map(e=>e.error).join('; '):String(g.control_overlap??'—'))])));
for(const l of DATA.loci){
 const card=el('section',undefined,'card');card.append(el('h2',`${l.genome_id} · ${l.locus_id}`));
 card.append(el('p',`${l.contig}:${l.start}–${l.end} (${l.strand}) · ${l.query_variants} variants / ${l.unique_sequences} unique sequences / ${l.support_species} labelled species`));
 const badge=el('span','','badge');card.append(badge);l.flags.forEach(f=>card.append(el('span',f,'badge bad')));
 card.append(el('p',`Proposed synteny window: ${l.contig}:${l.window_start}–${l.window_end}. Anchor boundaries: ${l.boundary_source}.`,'muted'));
 const scroll=el('div',undefined,'scroll');scroll.append(tracks(l));card.append(scroll);
 const toggle=el('button','Show tracks for every query variant','secondary');let fullTracks=false;
 toggle.addEventListener('click',()=>{fullTracks=!fullTracks;scroll.replaceChildren(tracks(l,fullTracks));toggle.textContent=fullTracks?'Show compact tracks':'Show tracks for every query variant'});card.append(toggle);
 card.append(el('p',`Representative BLAST chain: query coverage ${percent(l.coverage)}, identity ${percent(l.identity)}. Compact tracks show one chain per source gene and one predicted CDS model.`,'legend'));
 const detail=el('details');detail.append(el('summary','Reference support and alternative CDS models'));
 const support=[];for(const c of l.chains)for(const member of queries[c.query].members)support.push([
 member.query,member.gene_symbol||'',percent(c.coverage),percent(c.identity),c.evalue,c.hsps.length]);
 const wrap=el('div',undefined,'scroll');wrap.append(table(['Reference variant','Source gene','Query coverage','Identity','E-value','HSPs'],support));detail.append(wrap);
 card.append(detail);const row={review_id:DATA.review_id,genome_id:l.genome_id,locus_id:l.locus_id,decision:'',anchor_intervals:`${l.start}:${l.end}`,window_start:'',window_end:'',notes:''};
 const entry=decision(card,row,l);entry.badge=badge;
 for(const m of l.models){const p=el('p',`${m.model_id} · ${m.start}:${m.end} · coverage ${percent(m.coverage)} `),b=el('button','Use these anchor bounds','secondary');
 b.addEventListener('click',()=>{row.anchor_intervals=`${m.start}:${m.end}`;entry.anchor.value=row.anchor_intervals;persist()});p.append(b);detail.append(p)}
 if(l.refinement_error)detail.append(el('p',l.refinement_error,'error'));
 $('loci').append(card);
}
for(const g of DATA.genomes.filter(g=>!DATA.loci.some(l=>l.genome_id===g.genome_id))){const card=el('section',undefined,'card');card.append(el('h2',g.genome_id+' — '+g.status));
 card.append(el('p',g.error||'No retained locus. Record unresolved or exclude with a reason.'));
 decision(card,{review_id:DATA.review_id,genome_id:g.genome_id,locus_id:'',decision:'',anchor_intervals:'',window_start:'',window_end:'',notes:''},null);$('loci').append(card)}
function restore(rows){if(rows.length!==entries.length)throw Error('Decision row count differs from this checkpoint');
 const byKey=new Map(rows.map(r=>[r.genome_id+'\t'+r.locus_id,r]));if(byKey.size!==rows.length)throw Error('Duplicate decision rows');
 for(const e of entries){const r=byKey.get(e.row.genome_id+'\t'+e.row.locus_id);if(!r||r.review_id!==DATA.review_id)throw Error('Decision file belongs to another checkpoint')}
 for(const e of entries){Object.assign(e.row,byKey.get(e.row.genome_id+'\t'+e.row.locus_id));e.select.value=e.row.decision;if(e.anchor)e.anchor.value=e.row.anchor_intervals;e.notes.value=e.row.notes;
 const numbers=e.card.querySelectorAll('.decision label input');for(const n of numbers){const label=n.parentNode.firstChild.textContent;if(label==='Custom window start')n.value=e.row.window_start;if(label==='Custom window end')n.value=e.row.window_end}}
 persist();}
try{const saved=localStorage.getItem('gcf-stage1-'+DATA.review_id);if(saved)restore(JSON.parse(saved))}catch(_){}
function refresh(){const text=$('find').value.toLowerCase(),coverage=Number($('coverage').value),identity=Number($('identity').value),evalue=Number($('evalue').value);let visible=0,pass=0;
 for(const e of entries){const l=e.locus,passing=l&&l.chains.some(c=>!c.ambiguous_support&&c.coverage>=coverage&&c.identity>=identity&&c.evalue<=evalue);
 if(l){e.badge.textContent=passing?'Passes current display thresholds':'Below current display thresholds';e.badge.className='badge '+(passing?'good':'bad');if(passing)pass++}
 const match=(e.row.genome_id+' '+(l?l.locus_id+' '+l.contig:'')).toLowerCase().includes(text);
 const show=match&&(!$('passing').checked||passing);e.card.classList.toggle('hidden',!show);if(show)visible++}
 $('counts').textContent=`${DATA.loci.length} loci total · ${pass} pass current display thresholds · ${visible} review cards visible · ${entries.filter(e=>e.row.decision).length}/${entries.length} decisions entered`;
}
for(const id of ['find','coverage','identity','evalue','passing'])$(id).addEventListener('input',refresh);
$('loci').addEventListener('change',refresh);refresh();
$('export').addEventListener('click',()=>{const clean=x=>String(x??'').replace(/[\t\r\n]/g,' '),text=[columns.join('\t'),...entries.map(e=>columns.map(k=>clean(e.row[k])).join('\t'))].join('\n')+'\n';
 const url=URL.createObjectURL(new Blob([text],{type:'text/tab-separated-values;charset=utf-8'})),a=el('a');a.href=url;a.download='loci.decisions.tsv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)});
$('import').addEventListener('change',async event=>{try{const text=await event.target.files[0].text(),lines=text.replace(/[\r\n]+$/,'').split(/\r?\n/),headers=lines.shift().split('\t');
 restore(lines.map(line=>{const values=line.split('\t');return Object.fromEntries(headers.map((key,i)=>[key,values[i]??'']))}));refresh()}catch(error){alert(error.message)}});
</script></body></html>'''
