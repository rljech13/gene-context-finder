"""Keep gene-level GFF context, including features with no coding model."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from .sequence import open_text, sha256
from .early_synteny import attributes, datasets


def read_features(path, regions):
    by_contig=defaultdict(list)
    for r in regions:by_contig[r['contig']].append(r)
    features=[]
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('##FASTA'):break
            if line.startswith('#') or not line.strip():continue
            f=line.rstrip('\n').split('\t')
            if len(f)!=9:raise ValueError('Invalid GFF3 context row')
            if f[0] not in by_contig or not (f[2]=='gene' or f[2]=='pseudogene' or f[2].endswith('_gene')):continue
            a,b=int(f[3])-1,int(f[4]);attrs=attributes(f[8]);gid=attrs.get('ID')
            if not gid or not 0<=a<b:raise ValueError('Invalid annotated gene identity/coordinates')
            for r in by_contig[f[0]]:
                if a>=r['window_end'] or b<=r['window_start']:continue
                key=json.dumps([r['genome_id'],gid,f[0]],ensure_ascii=True)
                features.append(dict(feature_id=hashlib.sha256(key.encode()).hexdigest()[:24],
                    region_id=r['region_id'],genome_id=r['genome_id'],original_gene=gid,
                    name=attrs.get('gene') or attrs.get('Name') or gid,contig=f[0],start=a,end=b,
                    strand=f[6],kind=f[2],biotype=attrs.get('gene_biotype') or attrs.get('biotype') or f[2]))
    return features


def collect(data, manifest):
    """Enrich an existing checkpoint only from the exact original GFF release."""
    inputs={g['genome_id']:g for g in datasets(manifest)}
    receipts={g['genome_id']:g for g in data['provenance']['genomes']}
    features=[];sources=[]
    for gid in sorted({r['genome_id'] for r in data['regions']}):
        entry=inputs.get(gid,{})
        expected=receipts[gid].get('annotation_sha256')
        if not expected:
            sources.append(dict(genome_id=gid,status='no_original_annotation'));continue
        if not entry.get('gff3'):raise ValueError('Missing original GFF3 for '+gid)
        path=Path(entry['gff3']);before=path.stat();actual=sha256(path)
        if actual!=expected:raise ValueError('Annotation context differs from the original GFF3: '+gid)
        own=[r for r in data['regions'] if r['genome_id']==gid]
        found=read_features(path,own)
        if (path.stat().st_size,path.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):
            raise ValueError('GFF3 changed during context extraction: '+gid)
        features.extend(found);sources.append(dict(genome_id=gid,status='completed',gff3=str(path),sha256=actual,observations=len(found)))
    return dict(schema='stage2_annotation_context_v1',features=features,sources=sources,
        note='Gene-level annotation context, including noncoding genes. No CDS/model decisions or protein orthology are inferred.')
