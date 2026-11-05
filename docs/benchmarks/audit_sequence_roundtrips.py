"""Independent sequence/report checks using BioPython, never importing genecontext."""
import argparse, csv, hashlib, json, re
from collections import defaultdict, Counter
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Seq import Seq

p=argparse.ArgumentParser(); p.add_argument('comparison');p.add_argument('--source-code');p.add_argument('--out',required=True);a=p.parse_args()
r=Path(a.comparison).resolve();m=json.loads((r/'manifest.json').read_text());s=Path(m['source_results'])
def tab(path):
    with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))
def fa(path):return {x.id:str(x.seq) for x in SeqIO.parse(path,'fasta')}
def gff(path):
    models=defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line:continue
        f=line.split('\t')
        if f[2]!='CDS':continue
        attrs=dict(v.split('=',1) for v in f[8].split(';') if '=' in v)
        for parent in attrs['Parent'].split(','):models[parent].append((f[0],int(f[3])-1,int(f[4]),f[6],int(f[7]) if f[7]!='.' else 0))
    return models
def dna(parts,seqs,check_phase=True):
    parts=sorted(parts,key=lambda x:x[1],reverse=parts[0][3]=='-'); chunks=[]
    length=0
    for i,(contig,lo,hi,strand,phase) in enumerate(parts):
        piece=Seq(seqs[contig][lo:hi]);piece=piece.reverse_complement() if strand=='-' else piece
        if not i:piece=piece[phase:]
        elif check_phase:assert phase==(-length)%3,('phase',parts,i,length)
        chunks.append(str(piece));length+=len(piece)
    return ''.join(chunks).upper()
def translation(seq):return str(Seq(seq[:len(seq)//3*3]).translate()).removesuffix('*')
class Contig:
    def __init__(self,path,entry):self.path,self.entry=path,entry
    def __getitem__(self,region):
        length,offset,bases,width=self.entry;start=region.start or 0;end=length if region.stop is None else region.stop
        assert region.step in (None,1) and 0<=start<=end<=length
        if start==end:return ''
        lo=offset+(start//bases)*width+start%bases
        hi=offset+((end-1)//bases)*width+(end-1)%bases+1
        with open(self.path,'rb') as handle:handle.seek(lo);data=handle.read(hi-lo)
        result=data.replace(b'\n',b'').replace(b'\r',b'').decode();assert len(result)==end-start
        return result
class Genome:
    def __init__(self,path):
        self.path=Path(path);self.entries={}
        for line in Path(str(path)+'.fai').read_text().splitlines():
            f=line.split('\t');self.entries[f[0]]=tuple(map(int,f[1:5]))
    def __getitem__(self,contig):return Contig(self.path,self.entries[contig])
original=json.loads((s/'manifest.json').read_text()); genomes={}
for g in original['genomes']:
    local=s/'genomes'/g['genome_id']/'cache/genome.fa'
    genomes[g['genome_id']]=Genome(local if local.exists() else g['fasta'])
counts=Counter();windows={}
source_cds={q['id']:fa(s/'sequences'/q['id']/'cds.fna') for q in original['queries']}
source_aa={q['id']:fa(s/'sequences'/q['id']/'protein.faa') for q in original['queries']}
for row in tab(s/'models.tsv'):
    parts=[(row['contig'],lo,hi,row['strand'],phase or 0) for lo,hi,phase in json.loads(row['cds_exons'])]
    native=genomes[row['genome_id']];seq=dna(parts,native,check_phase=row['status']!='disrupted')
    expected=source_cds[row['query']][row['model_id']]
    if seq!=expected:
        extra=len(expected)-len(seq);assert 0<extra<=3
        terminal=min(p[1] for p in parts) if row['strand']=='-' else max(p[2] for p in parts)
        stop=native[row['contig']][terminal-extra:terminal] if row['strand']=='-' else native[row['contig']][terminal:terminal+extra]
        seq+=(str(Seq(stop).reverse_complement()) if row['strand']=='-' else stop).upper()
    assert seq==expected,(row['model_id'],'native CDS')
    # Source exports conservatively emit X for every non-ACGT codon, even
    # when BioPython can resolve e.g. CCN to Proline. Keep raw predictor
    # translation checks below unchanged; check this documented convention here.
    codons=[seq[i:i+3] for i in range(0,len(seq)-2,3)]
    unknown=sum(any(base not in 'ACGT' for base in codon) for codon in codons)
    source_translation=''.join('X' if any(base not in 'ACGT' for base in codon) else str(Seq(codon).translate()) for codon in codons).removesuffix('*')
    assert source_translation==source_aa[row['query']][row['model_id']],(row['model_id'],'protein translation')
    counts['source_ambiguous_codons_exported_as_X']+=unknown
    counts['source_models_with_ambiguous_codons']+=int(unknown>0)
    counts['search_models_native_roundtrip']+=1
for row in tab(r/'loci.tsv'):
    if row['comparison_error']:continue
    seq=genomes[row['genome_id']][row['contig']][int(row['window_start']):int(row['window_end'])]
    if row['strand']=='-':seq=str(Seq(seq).reverse_complement())
    path=r/row['window_fasta'];assert hashlib.sha256(path.read_bytes()).hexdigest()==row['window_sha256']
    assert fa(path)[row['locus_id']]==seq
    windows[row['locus_id']]=seq;counts['windows_native_roundtrip']+=1
for folder in sorted((r/'context_prediction').glob('tiberius_*')):
    models=gff(folder/'predictions.gff3');rawseq=fa(folder/'windows.fna')
    cds={k.split('|')[1]:v for k,v in fa(str(folder/'predictions.gff3')+'.cds.fna').items()}
    aa={k.split('|')[1]:v for k,v in fa(str(folder/'predictions.gff3')+'.proteins.faa').items()}
    assert set(models)==set(cds)==set(aa)
    for pid,parts in models.items():
        seq=dna(parts,rawseq);assert seq==cds[pid].upper();assert translation(seq)==aa[pid].removesuffix('*')
        counts['raw_models_roundtrip']+=1
cds=fa(r/'context_prediction_cds.fna');aa=fa(r/'context_prediction_proteins.faa')
for path in sorted((r/'context_predictions').glob('*.gff3')):
    for pid,parts in gff(path).items():
        pid=pid.removesuffix('.t1');seq=dna(parts,genomes[path.stem]);assert seq==cds[pid]
        assert translation(seq)==aa[pid];counts['native_models_roundtrip']+=1
for path in sorted((r/'queries').glob('*/cds.codon.aligned.fna')):
    proteins=fa(path.parent/'proteins.aligned.faa')
    for pid,seq in fa(path).items():
        codons=[seq[i:i+3] for i in range(0,len(seq),3)]
        assert all(c=='---' or '-' not in c for c in codons)
        aligned=''.join('-' if c=='---' else str(Seq(c).translate()) for c in codons)
        assert aligned==proteins[pid];counts['codon_alignments_roundtrip']+=1
class Links(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in {'href','src'} and value:
                url=urlsplit(value)
                if url.scheme or url.netloc or not url.path:continue
                target=(self.path.parent/unquote(url.path)).resolve()
                assert target.exists(),(str(self.path),value);counts['html_links']+=1
for path in [r/'report.html',*(r/'queries').glob('*/report.html')]:
    parser=Links();parser.path=path;parser.feed(path.read_text());counts['html_pages']+=1
for path in r.rglob('*.svg'):ET.parse(path);counts['svg_documents']+=1
if a.source_code:
    for name,digest in m['source_code_sha256'].items():
        assert hashlib.sha256((Path(a.source_code)/name).read_bytes()).hexdigest()==digest,name
        counts['source_modules']+=1
assert m['status']=='completed',m.get('errors')
assert all(b['status']=='completed' and b['cpu_execution']['visible_devices'] and not any('GPU' in d for d in b['cpu_execution']['visible_devices']) for b in m['context_prediction']['batches'])
counts['cpu_batches']=len(m['context_prediction']['batches'])
result=dict(status='passed',comparison=str(r),counts=dict(counts),errors=m['errors'],limits='Operation and sequence roundtrips, not gene-prediction accuracy or confirmed orthology.')
Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
