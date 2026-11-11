"""Read-only reuse of the STOP-2 genomic ruler and protein correspondences."""
import json
from pathlib import Path

from .stage2_report import source_pairs


def write_context(root, source, data, state, selected):
    root,source=Path(root),Path(source)
    features=data.get('annotation_features',[])
    supplement=source/'annotation_context.json'
    if supplement.exists():features=json.loads(supplement.read_text())['features']
    originals={r['region_id']:r for r in data['regions']};regions=[];annotation=[]
    for row in selected:
        r=dict(originals[row['source_region_id']],region_id=row['model_id'],source_region_id=row['source_region_id'],
            anchor_start=row['start'],anchor_end=row['end'],strand=row['strand'],copy_label=row['copy_label'],group=row['group'])
        regions.append(r)
        annotation.extend(dict(f,region_id=r['region_id']) for f in features if f['region_id']==r['source_region_id'])
    payload=dict(regions=regions,original_regions=data['regions'],models=data['models'],neighbor_pairs=data['neighbor_pairs'],
        source_neighbor_pairs=source_pairs(source,data),annotation_features=annotation,
        annotation_context_available=bool(features),decisions=state['models'],reviewed_pair_count=len(state.get('pairs',{})))
    (root/'context_evidence.json').write_text(json.dumps(payload,ensure_ascii=True,separators=(',',':'))+'\n')
    write_context_page(root,payload)


def write_context_page(root,payload):
    root=Path(root)
    template=Path(__file__).with_name('stage2_review.txt').read_text()
    css=template.split('<style>',1)[1].split('</style>',1)[0]
    helpers=Path(__file__).with_name('stage2_neighborhood.txt').read_text().split('function compare(){',1)[0]
    a=helpers.index('function openGeneModels(');b=helpers.index('// Coordinate view:',a)
    helpers=helpers[:a]+helpers[b:]
    cluster=Path(__file__).with_name('stage2_cluster.txt').read_text()
    code=r'''
'use strict';const DATA=__DATA__,$=id=>document.getElementById(id),regions=new Map([...DATA.original_regions,...DATA.regions].map(r=>[r.region_id,r])),models=new Map(DATA.models.map(m=>[m.model_id,m]));
const wantedGroup=new URLSearchParams(location.search).get('group');if(wantedGroup)DATA.regions=DATA.regions.filter(r=>r.group===wantedGroup);
function el(t,text,cls){const n=document.createElement(t);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n}
function svgNode(t,a,text){const n=document.createElementNS('http://www.w3.org/2000/svg',t);for(const[k,v]of Object.entries(a))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n}
function current(){return regions.get($('region').value)}
function own(r){return DATA.models.filter(m=>m.region_id===r.source_region_id).map(m=>({...m,overlaps_anchor:m.strand===r.strand&&m.start<r.anchor_end&&m.end>r.anchor_start}))}
function groupModels(r){const out=new Map();for(const m of own(r)){if(!out.has(m.gene_group))out.set(m.gene_group,[]);out.get(m.gene_group).push(m)}return [...out.values()].sort((a,b)=>Math.min(...a.map(m=>m.start))-Math.min(...b.map(m=>m.start)))}
function groupLabel(g){const names=[];for(const m of g){for(const n of [m.identification.proposed_name,m.annotation_name,...m.identification.near_best.map(h=>h.gene_symbol)].filter(Boolean))if(!names.some(x=>x.toUpperCase()===n.toUpperCase()))names.push(n)}return names.join(' / ')||'Unassigned gene'}
__HELPERS__
__CLUSTER__
function openGeneModels(group){const r=current(),g=own(r).filter(m=>m.gene_group===group),box=$('modelInfo');box.replaceChildren(el('h2',groupLabel(g)));for(const m of g){const c=DATA.decisions[m.model_id]||{};box.append(el('p',`${m.source} · ${m.original_model} · ${m.protein_length} aa · ${m.contig}:${m.start}–${m.end} (${m.strand}) · ${c.decision==='inspect'?'reviewed '+c.role:'context evidence, not selected'}`))}box.hidden=false;box.scrollIntoView({behavior:'smooth',block:'nearest'})}
function compare(){if(!current())return;const rows=comparisonRows(),edges=rows.map((r,i)=>i?correspondenceEdges(rows[0],r):[]);diagramFilterNotice(rows);$('links').replaceChildren(neighborhoodDiagram(rows,edges));$('modelInfo').hidden=true}
for(const r of DATA.regions){const o=el('option',regionLabel(r));o.value=r.region_id;$('region').append(o)}
function refreshOther(){const old=$('other').value;$('other').replaceChildren();const all=el('option','All selected copies');all.value='all';$('other').append(all);for(const r of DATA.regions.filter(r=>r.region_id!==current()?.region_id)){const o=el('option',regionLabel(r));o.value=r.region_id;$('other').append(o)}if([...$('other').options].some(o=>o.value===old))$('other').value=old;compare()}
$('region').addEventListener('change',refreshOther);for(const id of ['other','scale','neighborsShown','genomeZoom'])$(id).addEventListener('change',compare);
const wanted=new URLSearchParams(location.search).get('model');if(DATA.regions.some(r=>r.region_id===wanted))$('region').value=wanted;
$('reviewedStatus').textContent=`${Object.values(DATA.decisions).filter(c=>c.decision==='inspect'&&c.role==='neighbor').length} manually selected neighboring model(s) · ${DATA.reviewed_pair_count} reviewed cross-locus correspondence(s). Other arrows and ribbons remain descriptive evidence.`;
refreshOther();
'''
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    code=code.replace('__DATA__',text).replace('__HELPERS__',helpers).replace('__CLUSTER__',cluster)
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Available gene context — read-only</title><style>'''+css+'''</style><main>
<h1>Available gene context</h1><p><a href="report.html">Back to the reference sheet</a></p>
<p class="muted">Read-only annotation + Tiberius evidence from STOP 2. Arrows use a shared kb ruler; zero is the selected CDS-model midpoint. Ribbons are saved protein-match proposals, not confirmed orthology or DNA alignments.</p>
<p id="reviewedStatus" class="muted"></p><section>
<div class="controls"><label>Selected copy <select id="region"></select></label><label>Compare with <select id="other"></select></label></div>
<div class="controls display-controls"><label>Show on diagram <select id="modelSource"><option value="all">Annotation + Tiberius</option><option value="annotation">Original annotation</option><option value="tiberius">Tiberius predictions</option></select></label><label>View <select id="scale"><option value="distance">Gene clusters (kb)</option><option value="models">CDS model details</option><option value="order">Gene order (schematic)</option></select></label><label>Zoom <select id="genomeZoom"><option value="1">1×</option><option value="2">2×</option><option value="4">4×</option></select></label><label>Neighbors per side <select id="neighborsShown"><option value="5">Up to 5</option><option value="10" selected>Up to 10</option></select></label></div>
<div id="diagramStatus" class="filter-notice" hidden></div><div id="links" class="scroll"></div><p class="muted">A+T: identical CDS alternatives · A/T: differing alternatives from both sources · A: annotation · T: Tiberius · RNA: noncoding GFF context. Click a gene in the first row to inspect sources. No decisions can be changed here.</p>
</section><section id="modelInfo" hidden></section></main><script>'''+code+'</script></html>'
    (root/'context.html').write_text(page)
