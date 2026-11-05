"""Independent full Swiss-Prot comparison for candidates and window predictions.
Uses BLAST+ directly, never imports genecontext. Top hits are evidence, not orthology.
"""
from pathlib import Path
import argparse,csv,gzip,json,re,subprocess
from collections import defaultdict
from html import escape
p=argparse.ArgumentParser();p.add_argument('--comparison');p.add_argument('--reference',required=True);p.add_argument('--work',required=True);p.add_argument('--threads',default='16');p.add_argument('--prepare',action='store_true');a=p.parse_args()
w=Path(a.work);w.mkdir(parents=True,exist_ok=True)
def fasta(path):
 records={};key=None;chunks=[]
 opener=gzip.open if str(path).endswith('.gz') else open
 with opener(path,'rt') as h:
  for line in h:
   if line.startswith('>'):
    if key is not None:records[key]=''.join(chunks)
    key=line[1:].split()[0];chunks=[]
   else:chunks.append(line.strip())
 if key is not None:records[key]=''.join(chunks)
 return records
def table(path):
 with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))
def union(parts):
 length=0;end=-1
 for lo,hi in sorted(parts):
  if hi>end:length+=hi-max(lo,end);end=hi
 return length
library=w/'swissprot.faa';meta_path=w/'swissprot.headers.json';db=w/'swissprot'
if not library.exists():
 meta={};opener=gzip.open if a.reference.endswith('.gz') else open
 with opener(a.reference,'rt') as inp,open(library,'w') as out:
  for line in inp:
   out.write(line)
   if line.startswith('>'):
    header=line[1:].strip();identifier=header.split()[0]
    gene=re.search(r'\bGN=(\S+)',header);organism=re.search(r'\bOS=(.+?)(?: OX=| GN=| PE=| SV=|$)',header)
    meta[identifier]={'description':header,'gene':gene.group(1) if gene else '', 'organism':organism.group(1) if organism else ''}
 meta_path.write_text(json.dumps(meta))
if not Path(str(db)+'.pin').exists():
 with open(w/'makeblastdb.log','w') as log:subprocess.run(['makeblastdb','-in',str(library),'-dbtype','prot','-out',str(db)],stdout=log,stderr=subprocess.STDOUT,check=True)
if a.prepare:print('Swiss-Prot audit database ready');raise SystemExit(0)
r=Path(a.comparison).resolve();out=r/'validation';out.mkdir(exist_ok=True)
m=json.loads((r/'manifest.json').read_text());source=Path(m['source_results']);original=json.loads((source/'manifest.json').read_text())
meta=json.loads(meta_path.read_text());sequences={};info={}
source_aa={q['id']:fasta(source/'sequences'/q['id']/'protein.faa') for q in original['queries']}
for row in table(r/'candidates.tsv'):
 if row['selected']!='True' and not (row['representative']=='True' and int(row['locus_rank'])<=3):continue
 qid=row['query'];mid=row['model_id']
 if mid not in source_aa[qid]:continue
 key='C'+str(len(info)+1);sequences[key]=source_aa[qid][mid]
 info[key]={'kind':'selected_target' if row['selected']=='True' else 'alternative_locus','genome_id':row['genome_id'],'group':row['group'],'source_gene':row['source_gene'],'model_id':mid,'protein_aa':len(sequences[key])}
raw=fasta(r/'context_prediction_proteins.faa')
for row in table(r/'context_predictions.tsv'):
 pid=row['protein_id'];key='P'+str(len(info)+1);sequences[key]=raw[pid]
 info[key]={'kind':'Tiberius_target_overlap' if row['target_overlap']=='True' else 'Tiberius_neighbor','genome_id':row['genome_id'],'group':'','source_gene':row.get('name',''),'model_id':pid,'protein_aa':len(sequences[key])}
queries=out/'identity_audit_queries.faa'
with open(queries,'w') as h:
 for key,seq in sequences.items():h.write('>'+key+'\n'+seq+'\n')
hits=out/'identity_audit.raw.tsv'
cmd=['blastp','-query',str(queries),'-db',str(db),'-num_threads',a.threads,'-evalue','1e-5','-max_target_seqs','10','-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen']
with open(hits,'w') as h,open(out/'identity_audit.blastp.log','w') as log:subprocess.run(cmd,stdout=h,stderr=log,check=True)
by_hit=defaultdict(list)
for line in hits.read_text().splitlines():
 f=line.split('\t');by_hit[f[0],f[1]].append(f)
options=defaultdict(list)
for (qid,sid),hsps in by_hit.items():
 qa=[];sa=[];chosen=[]
 for h in sorted(hsps,key=lambda h:-float(h[9])):
  q=(int(h[4])-1,int(h[5]));s=(int(h[6])-1,int(h[7]))
  if any(max(q[0],x)<min(q[1],y) for x,y in qa) or any(max(s[0],x)<min(s[1],y) for x,y in sa):continue
  qa.append(q);sa.append(s);chosen.append(h)
 n=sum(int(h[3]) for h in chosen);ident=sum(float(h[2])*int(h[3]) for h in chosen)/n
 options[qid].append({'swissprot_id':sid,**meta.get(sid,{}),'bitscore':sum(float(h[9]) for h in chosen),'evalue':min(float(h[8]) for h in chosen),'identity_percent':ident,'query_coverage':union(qa)/len(sequences[qid]),'reference_coverage':union(sa)/int(chosen[0][11])})
rows=[]
for qid,item in info.items():
 opts=sorted(options[qid],key=lambda d:-d['bitscore'])[:5]
 if not opts:rows.append({'query_id':qid,**item,'rank':0,'status':'no_significant_swissprot_hit'})
 for rank,hit in enumerate(opts,1):rows.append({'query_id':qid,**item,'rank':rank,'status':'significant_local_homology',**hit})
cols=['query_id','kind','genome_id','group','source_gene','model_id','protein_aa','rank','status','swissprot_id','gene','organism','bitscore','evalue','identity_percent','query_coverage','reference_coverage','description']
with open(out/'identity_audit.tsv','w') as h:
 writer=csv.DictWriter(h,cols,delimiter='\t');writer.writeheader();writer.writerows(rows)
page=['<!doctype html><meta charset="utf-8"><title>Independent protein identity audit</title><style>body{font:14px system-ui;margin:30px}table{border-collapse:collapse}td,th{border:1px solid #ccd;padding:7px}</style><h1>Independent protein identity audit</h1><p>Top Swiss-Prot BLASTP matches for selected targets, alternative loci and every Tiberius prediction. Local/domain matches are not individual-copy orthology or validated gene names.</p><p><a href="report.html">Window prediction audit</a> · <a href="identity_audit.tsv">All top-five matches</a></p><table><tr>']
shown=['kind','genome_id','source_gene','gene','organism','identity_percent','query_coverage','reference_coverage','evalue']
page += [''.join('<th>'+escape(k)+'</th>' for k in shown)+'</tr>']
for row in rows:
 if row['rank'] not in (0,1):continue
 page.append('<tr>'+''.join('<td>'+escape(f'{row.get(k,""):.3f}' if isinstance(row.get(k),float) and k!='evalue' else str(row.get(k,'')))+'</td>' for k in shown)+'</tr>')
page.append('</table>');(out/'identity_audit.html').write_text('\n'.join(page))
(out/'identity_audit.json').write_text(json.dumps({'queries':len(info),'command':cmd,'reference':a.reference,'results':len(rows),'limits':'Swiss-Prot best hits do not independently establish orthology; inspect coverage, competing loci and annotation.'},indent=2))
print('Audited',len(info),'candidate/predicted proteins')
