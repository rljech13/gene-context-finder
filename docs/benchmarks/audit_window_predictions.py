"""Independent CDS-overlap audit of window predictions against supplied GFF3.
No genecontext imports. Agreement is not a biological accuracy estimate.
"""
import argparse, csv, gzip, json
from pathlib import Path
from collections import defaultdict, Counter
from html import escape
from urllib.parse import unquote

p=argparse.ArgumentParser();p.add_argument('comparison');p.add_argument('--panel',required=True);a=p.parse_args()
r=Path(a.comparison).resolve();out=r/'validation';out.mkdir(exist_ok=True)
def table(f):
 with open(f) as h:return list(csv.DictReader(h,delimiter='\t'))
def write_table(f,rs,cols):
 with open(f,'w') as h:
  w=csv.DictWriter(h,cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rs)
def merged(parts):
 result=[]
 for lo,hi in sorted(parts):
  if result and lo<=result[-1][1]:result[-1]=(result[-1][0],max(result[-1][1],hi))
  else:result.append((lo,hi))
 return result
def length(parts):return sum(b-a for a,b in merged(parts))
def overlap(a,b):return length([(max(x,u),min(y,v)) for x,y in a for u,v in b if max(x,u)<min(y,v)])
def parse_gff(path,windows=None):
 models={};genes={}
 with open(path,'rb') as h:compressed=h.read(2)==b'\x1f\x8b'
 opener=gzip.open if compressed else open
 with opener(path,'rt') as h:
  for line in h:
   if not line.strip() or line.startswith('#'):continue
   f=line.rstrip().split('\t')
   if len(f)!=9:continue
   contig,lo,hi,strand=f[0],int(f[3])-1,int(f[4]),f[6]
   if windows is not None and not any(contig==w['contig'] and lo<int(w['window_end']) and hi>int(w['window_start']) for w in windows):continue
   attrs={k:unquote(v) for k,v in (piece.split('=',1) for piece in f[8].split(';') if '=' in piece)}
   if f[2]=='gene':genes[attrs.get('ID','')]={'name':attrs.get('Name',attrs.get('gene','unnamed')),'start':lo,'end':hi}
   if f[2] in ('mRNA','transcript'):
    mid=attrs.get('ID','')
    model=models.setdefault(mid,{'id':mid,'parts':[]})
    model.update(gene=attrs.get('Parent',mid),contig=contig,strand=strand,start=lo,end=hi)
   if f[2]=='CDS':
    for mid in attrs.get('Parent','').split(','):
     if not mid:continue
     model=models.setdefault(mid,{'id':mid,'parts':[]})
     model['parts'].append((lo,hi));model.update(contig=contig,strand=strand)
 for model in models.values():
  if model['parts']:
   model.setdefault('start',min(x for x,y in model['parts']));model.setdefault('end',max(y for x,y in model['parts']))
   model.setdefault('gene',model['id']);model['name']=genes.get(model['gene'],{}).get('name',model['gene'])
 return [m for m in models.values() if m['parts']],genes

def compare_model(pred,ref):
 if pred['contig']!=ref['contig'] or pred['strand']!=ref['strand']:return (0,0,0)
 common=overlap(pred['parts'],ref['parts']);lp=length(pred['parts']);lr=length(ref['parts'])
 return (2*common/(lp+lr),common/lr,common/lp)
manifest=json.loads((r/'manifest.json').read_text());source=Path(manifest['source_results']);original=json.loads((source/'manifest.json').read_text())
panel=table(a.panel);inputs={g['genome_id']:g for g in original['genomes']}
loci=table(r/'loci.tsv');summary=table(r/'summary.tsv');predictions=table(r/'context_predictions.tsv');neighbors=table(r/'neighbors.tsv')
by_genome=defaultdict(list)
for row in loci:
 if not row.get('comparison_error'):by_genome[row['genome_id']].append(row)
by_pid={v['protein_id']:v for v in predictions}
records=[];model_records=[];figures=[]
for gid,windows in by_genome.items():
 entry=inputs[gid];has_annotation=bool(entry.get('gff3'))
 refs,_=parse_gff(entry['gff3'],windows) if has_annotation else ([],{})
 native=r/'context_predictions'/(gid+'.gff3');preds,_=parse_gff(native) if native.exists() else ([],{})
 for w in windows:
  lid=w['locus_id'];lo=int(w['window_start']);hi=int(w['window_end']);contig=w['contig']
  ref=[m for m in refs if m['contig']==contig and m['start']<hi and m['end']>lo]
  ps=[m for m in preds if m['id'].removesuffix('.t1').startswith(lid+'.')]
  for m in ps:
   meta=by_pid[m['id'].removesuffix('.t1')];m['name']=meta.get('name') or ('target-overlapping prediction' if meta['target_overlap']=='True' else m['id'])
  full=defaultdict(list)
  for m in ref:
   if m['start']>=lo and m['end']<=hi:full[m['gene']].append(m)
  matched=0
  for gene,alternatives in full.items():
   scores=[(*compare_model(pm,rm),pm['id'],rm['id']) for pm in ps for rm in alternatives]
   best=max(scores,default=(0,0,0,'',''))
   passed=best[1]>=0.8 and best[2]>=0.8;matched+=passed
   model_records.append(dict(locus_id=lid,genome_id=gid,annotation_gene=gene,annotation_name=alternatives[0]['name'],reference_model=best[4],predicted_model=best[3],cds_dice=best[0],reference_cds_recall=best[1],prediction_cds_precision=best[2],agreement_80pct=passed))
  target_refs=[m for m in ref if m['id']==w.get('source_model')]
  ts=max([(*compare_model(pm,rm),pm['id']) for pm in ps for rm in target_refs],default=(0,0,0,''))
  record=dict(locus_id=lid,group=w['group'],genome_id=gid,target=w['target_name'],annotation_available=has_annotation,
   tiberius_status=w['context_status'],model=w['context_predictor_model'],predicted_models=len(ps),predicted_neighbors=sum(v['locus_id']==lid for v in neighbors),
   predicted_target_overlap=sum(v['locus_id']==lid and v['target_overlap']=='True' for v in predictions),
   annotated_genes_in_window=len({m['gene'] for m in ref}),fully_contained_annotated_coding_genes=len(full),matched_annotated_genes_80pct=matched,
   annotation_gene_agreement=matched/len(full) if full else '',target_reference_cds_recall=ts[1] if target_refs else '',target_prediction_cds_precision=ts[2] if target_refs else '',
   figure=lid+'.svg')
  records.append(record)
  # Native-coordinate CDS tracks; one longest reference isoform per gene is drawn.
  # Agreement above evaluates every overlapping isoform, not just the drawn one.
  by_gene=defaultdict(list)
  for m in ref:by_gene[m['gene']].append(m)
  shown=[max(ms,key=lambda m:length(m['parts'])) for ms in by_gene.values()]
  tracks=[('Supplied GFF3 — representative coding isoform',shown,'#2563eb'),('Tiberius — DNA-only predictions',ps,'#059669')]
  shapes=['<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="HEIGHT" viewBox="0 0 1600 HEIGHT">','<rect width="100%" height="100%" fill="white"/>','<style>text{font:12px system-ui;fill:#172033}</style>',f'<text x="20" y="26" font-size="18">{escape(gid)} | {escape(w["target_name"])} | {lid}</text>',f'<text x="20" y="50">Native window: {escape(contig)}:{lo}–{hi} (0-based, half-open); target strand {w["strand"]}; boxes = CDS exons, lines = introns</text>']
  scale=1460/(hi-lo);xx=lambda pos:100+(max(lo,min(hi,pos))-lo)*scale;y=90
  tx1,tx2=xx(int(w['target_start'])),xx(int(w['target_end']))
  # A narrow target guide, shared across both evidence tracks.
  for title,models,color in tracks:
   shapes.append(f'<text x="20" y="{y}">{escape(title)}</text>');y+=28
   if not models:
    shapes.append(f'<text x="100" y="{y}">No models in this track.</text>');y+=28
   for model in sorted(models,key=lambda m:(m['start'],m['id'])):
    c='#dc2626' if (model.get('gene')==w.get('source_gene') or by_pid.get(model['id'].removesuffix('.t1'),{}).get('target_overlap')=='True') else color
    first,last=xx(min(x for x,z in model['parts'])),xx(max(z for x,z in model['parts']))
    shapes.append(f'<path d="M{first:.2f} {y}H{last:.2f}" stroke="{c}"/>')
    for x,z in model['parts']:
     if x<hi and z>lo:shapes.append(f'<rect x="{xx(x):.2f}" y="{y-5}" width="{max(1,xx(z)-xx(x)):.2f}" height="10" fill="{c}"/>')
    direction='→' if model['strand']=='+' else '←'
    label=f'{model["name"]} {direction} | {model["id"]}'
    shapes.append(f'<text x="100" y="{y+18}" style="font-size:10px">{escape(label)}</text>');y+=42
   y+=18
  for tick in range(6):
   pos=lo+(hi-lo)*tick/5
   shapes.append(f'<text x="{xx(pos):.2f}" y="{y}" text-anchor="middle">{int(pos):,}</text>')
  shapes.append(f'<text x="20" y="{y+28}">Red = selected target / target-overlapping prediction. Predictions and names remain hypotheses. Native GFF3 coordinates are 1-based inclusive.</text>')
  shapes.append('</svg>');height=y+50
  (out/record['figure']).write_text('\n'.join(shapes).replace('HEIGHT',str(height)))

cols=list(records[0]) if records else ['locus_id'];write_table(out/'window_prediction_audit.tsv',records,cols)
write_table(out/'annotation_gene_agreement.tsv',model_records,list(model_records[0]) if model_records else ['locus_id'])
report=['<!doctype html><meta charset="utf-8"><title>Pipeline validation</title><style>body{font:15px system-ui;max-width:1700px;margin:30px auto;padding:20px}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #ccd;padding:8px;text-align:left}img{width:100%}small{color:#475569}</style><h1>Pipeline validation</h1>',
 '<p><a href="../report.html">Comparative report and alignments</a> · <a href="window_prediction_audit.tsv">Window prediction audit</a> · <a href="annotation_gene_agreement.tsv">CDS agreement per annotation gene</a> · <a href="../context_prediction.json">Tiberius runtime and model provenance</a></p>',
 '<p>The supplied annotation and DNA-only Tiberius predictions are compared independently. A reference gene agrees when at least one predicted model covers ≥80% of its annotated CDS and ≥80% of that prediction overlaps the same reference CDS on the same strand. All overlapping annotated isoforms are considered. Only genes fully contained within the window enter this fraction. Annotation itself is not ground truth; this measures agreement, not biological sensitivity or specificity.</p>',
 '<table><tr><th>Species / assembly level</th><th>Query</th><th>Selected candidate</th><th>Search status</th><th>Tiberius</th><th>Annotation agreement</th><th>Evidence</th></tr>']
for sp in panel:
 gid=sp['genome_id']
 for group in manifest['groups']:
  obs=next((v for v in summary if v['genome_id']==gid and v['group']==group['id']),{})
  rec=next((v for v in records if v['genome_id']==gid and v['group']==group['id']),{})
  link=f'<a href="{rec["figure"]}">Exon tracks</a> · <a href="../context_predictions/{escape(gid)}.gff3">Native prediction GFF3</a>' if rec and (r/'context_predictions'/(gid+'.gff3')).exists() else 'No accessible prediction track'
  agree=f'{rec["matched_annotated_genes_80pct"]}/{rec["fully_contained_annotated_coding_genes"]} fully contained genes' if rec.get('annotation_available') else 'No supplied GFF3'
  pred=f'{rec["predicted_models"]} models; {rec["predicted_neighbors"]} neighbours; {rec["predicted_target_overlap"]} target-overlapping' if rec else 'No window / inspect search status'
  report.append('<tr>'+''.join('<td>'+v+'</td>' for v in [escape(sp['species'])+'<br><small>'+escape(sp['accession']+'; '+sp['assembly_level'])+'</small>',escape(group['name']),escape(obs.get('target_name','')),escape(obs.get('status','not_processed')),escape(pred),escape(agree),link])+'</tr>')
report+=['</table><h2>Window details</h2>']
for rec in records:
 report += [f'<h3>{escape(rec["genome_id"])} — {escape(rec["target"])}</h3>',f'<p>Tiberius model: {rec["model"]}; status: {rec["tiberius_status"]}; target CDS recall/precision against supplied annotation: {rec["target_reference_cds_recall"]} / {rec["target_prediction_cds_precision"]}.</p>',f'<img src="{rec["figure"]}" alt="Supplied and predicted CDS exon tracks">']
(out/'report.html').write_text('\n'.join(report))
result={'status':'completed','input_genomes':len(inputs),'requested_genomes':len(panel),'windows':len(records),'prediction_models':len(predictions),'fully_contained_reference_genes':len(model_records),'matching_reference_genes_80pct':sum(x['agreement_80pct'] for x in model_records),'limits':'CDS agreement with supplied annotation; not ground truth, orthology, or gene absence.'}
(out/'audit_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
