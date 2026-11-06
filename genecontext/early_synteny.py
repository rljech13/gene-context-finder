"""Automatic, cached annotation-supported local context evidence for stage 1.

Reference origin is resolved by versioned identifiers AND exact protein sequence,
then by exact sequence when identifiers are unavailable. Similarity is never used
to invent a reference origin. All correspondences remain tentative homology.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
import time
from urllib.parse import unquote
from uuid import uuid4

from .sequence import Model, coding_sequence, open_text, sha256, translate, write_fasta, write_tsv
from .tools import Genome, Runner, tool_versions

SCHEMA = 'early_synteny_v1'
ACCESSION = re.compile(r'(?<![A-Za-z0-9])(?:[NXY][PMR]_\d+\.\d+)(?![A-Za-z0-9.])')
AA = set('ACDEFGHIKLMNPQRSTVWYXBZJUO')


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def settings(args):
    return dict(neighbors=getattr(args, 'synteny_neighbors', 5),
                max_neighbors=getattr(args, 'synteny_max_neighbors', 20),
                max_distance=getattr(args, 'synteny_max_distance', 1000000),
                min_coverage=.5, min_identity=.3, competing_ratio=.95, evalue=1e-5)


def datasets(manifest):
    """Resolve all dataset paths relative to their own manifest, never cwd."""
    manifest = Path(manifest).resolve()
    with manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if 'genome_id' not in (reader.fieldnames or []):
            raise ValueError('Dataset TSV requires genome_id')
        rows = []
        for source in reader:
            if not any(source.values()):
                continue
            row = dict(source)
            if not row.get('genome_id'):
                raise ValueError('Dataset requires a nonempty genome_id')
            for key in ['fasta', 'gff3', 'protein', 'cds', 'annotation_provenance']:
                row[key] = str((manifest.parent/row[key]).resolve()) if row.get(key) else ''
            # Only exact NCBI sibling names are discovered; no recursive folder
            # guessing and no mixing annotation releases from other directories.
            if not row['gff3'] and row['fasta']:
                fasta = Path(row['fasta'])
                prefix = re.sub(r'_genomic\.(?:fna|fa|fasta)(?:\.gz)?$', '', fasta.name)
                if prefix != fasta.name:
                    choices = [fasta.parent/(prefix+'_genomic.gff'+suffix) for suffix in ['', '.gz']]
                    row['gff3'] = next((str(p) for p in choices if p.is_file()), '')
            if row['gff3']:
                gff = Path(row['gff3'])
                prefix = re.sub(r'_genomic\.gff(?:\.gz)?$', '', gff.name)
                if prefix != gff.name:
                    if not row['protein']:
                        row['protein'] = next((str(p) for p in [gff.parent/(prefix+'_protein.faa'+s) for s in ['', '.gz']] if p.is_file()), '')
                    if not row['cds']:
                        row['cds'] = next((str(p) for p in [gff.parent/(prefix+'_cds_from_genomic.fna'+s) for s in ['', '.gz']] if p.is_file()), '')
            rows.append(row)
    ids = [r['genome_id'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate dataset genome_id')
    return rows


def attributes(raw):
    return {unquote(k): unquote(v) for item in raw.split(';') if '=' in item
            for k, v in [item.split('=', 1)]}


def parse_annotation(path):
    genes, meta, models, protein_ids, regions = {}, {}, {}, defaultdict(set), {}
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('##sequence-region '):
                _, contig, begin, end = line.split()
                regions[contig] = (int(begin)-1, int(end))
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Invalid GFF3 row {number}')
            contig, _, kind, start, end, _, strand, phase, raw = fields
            start, end = int(start)-1, int(end)
            if not 0 <= start < end:
                raise ValueError('Invalid GFF3 coordinates')
            attrs = attributes(raw); fid = attrs.get('ID', '')
            if fid:
                if fid in meta and kind not in {'CDS', 'exon', 'stop_codon'}:
                    raise ValueError('Repeated annotation ID: '+fid)
                meta[fid] = dict(attrs, _contig=contig, _strand=strand)
            if kind == 'gene':
                if not fid:
                    raise ValueError('Gene without ID')
                genes[fid] = dict(gene_id=fid, name=attrs.get('Name', fid), contig=contig,
                    start=start, end=end, strand=strand, biotype=attrs.get('gene_biotype', ''),
                    ncbi_gene_id=next((x[7:] for x in attrs.get('Dbxref', '').split(',') if x.startswith('GeneID:')), ''),
                    models=[])
            if kind not in {'CDS', 'stop_codon'}:
                continue
            parents = attrs.get('Parent', '').split(',')
            if not parents[0] or strand not in {'+', '-'}:
                raise ValueError('Coding feature requires Parent and strand')
            value = None if phase == '.' else int(phase)
            if value not in {None, 0, 1, 2}:
                raise ValueError('Invalid GFF3 phase')
            for parent in parents:
                model = models.setdefault(parent, Model(parent, contig, strand))
                if (model.contig, model.strand) != (contig, strand):
                    raise ValueError('Transcript spans contigs or strands: '+parent)
                if kind == 'CDS':
                    model.cds.append((start, end, value))
                    if attrs.get('protein_id'):
                        protein_ids[parent].add(attrs['protein_id'])
                else:
                    model.stops.append((start, end))
    retained = {}
    for mid, model in models.items():
        if not model.cds:
            continue
        attrs = meta.get(mid, {})
        parent = mid if mid in genes else attrs.get('Parent', '')
        if ',' in parent:
            raise ValueError('Transcript assigned to several genes: '+mid)
        if parent not in genes:
            # Do not silently invent gene order from disconnected CDS features.
            continue
        gene = genes[parent]
        if model.contig != gene['contig'] or model.strand != gene['strand']:
            raise ValueError('Transcript/gene coordinate identity differs: '+mid)
        if any(a < gene['start'] or b > gene['end'] for a, b, _ in model.cds):
            raise ValueError('CDS outside parent gene: '+mid)
        model.gene = parent; model.attrs = attrs
        gene['models'].append(mid)
        retained[mid] = dict(model=model,
            protein_ids=sorted(protein_ids[mid] | ({attrs['protein_id']} if attrs.get('protein_id') else set())),
            aliases=sorted(set([mid, mid.removeprefix('rna-'), attrs.get('transcript_id', ''), attrs.get('protein_id', '')]) | protein_ids[mid]))
    for gene in genes.values():
        gene['coding'] = bool(gene['models'])
    return genes, retained, regions


def fasta_records(path):
    header, pieces = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(pieces).upper().removesuffix('*')
                header, pieces = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError('Sequence before FASTA header')
                pieces.append(''.join(line.split()))
        if header is not None:
            yield header, ''.join(pieces).upper().removesuffix('*')


@contextmanager
def cache_lock(path):
    """OS-released advisory lock; a killed worker cannot leave a stale lock."""
    import fcntl
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def verified_cache(path, identity):
    if not (path/'cache.json').is_file():
        return False
    record = json.loads((path/'cache.json').read_text())
    if record.get('identity') != identity:
        raise ValueError('Cache input identity differs')
    for name, digest in record['outputs'].items():
        f = (path/name).resolve()
        if path.resolve() not in f.parents or not f.is_file() or sha256(f) != digest:
            raise ValueError('Cache evidence changed: '+name)
    return True


def seal_cache(path, identity):
    save(path/'cache.json', dict(identity=identity,
        outputs={str(p.relative_to(path)): sha256(p) for p in sorted(path.rglob('*')) if p.is_file() and p.name != 'cache.json'}))


def build_catalog(row, cache_base, runner, known_fasta_digest=None):
    t0 = time.perf_counter()
    gff = Path(row['gff3'])
    use_genome=not (row.get('protein') or row.get('cds'))
    fasta=Path(row['fasta']) if row.get('fasta') else None
    if use_genome and (fasta is None or not fasta.is_file()):
        raise ValueError('Need supplied protein/CDS FASTA or a genomic FASTA to derive products')
    input_paths=[gff]+([fasta] if use_genome else [])+[Path(row[k]) for k in ['protein','cds','annotation_provenance'] if row.get(k)]
    initial_stats={str(p):(p.stat().st_size,p.stat().st_mtime_ns) for p in input_paths}
    inputs = {'gff3': dict(path=str(gff), sha256=sha256(gff))}
    if use_genome:
        inputs['fasta']=dict(path=str(fasta),sha256=known_fasta_digest or sha256(fasta))
    if row.get('protein'):
        inputs['protein'] = dict(path=row['protein'], sha256=sha256(row['protein']))
    if row.get('cds'):
        inputs['cds'] = dict(path=row['cds'], sha256=sha256(row['cds']))
    if row.get('annotation_provenance'):
        provenance_path=Path(row['annotation_provenance'])
        provenance=json.loads(provenance_path.read_text())
        if provenance.get('genome_id')!=row['genome_id'] or provenance.get('snapshot_sha256')!=inputs['gff3']['sha256']:
            raise ValueError('Annotation snapshot provenance does not match this dataset')
        if provenance.get('source_sha256')!=provenance.get('expected_sha256') or not provenance.get('source_sha256'):
            raise ValueError('Original annotation fingerprint is not verified')
        inputs['annotation_provenance']=dict(path=str(provenance_path),sha256=sha256(provenance_path))
        inputs['gff3']['derivation']=provenance
    identity = dict(schema=SCHEMA, implementation_sha256=sha256(Path(__file__)), inputs={k:v['sha256'] for k,v in inputs.items()})
    key = token(identity); destination = cache_base/'catalogs'/key
    with cache_lock(destination.with_suffix('.lock')):
        reused = verified_cache(destination, identity)
        if not reused:
            tmp = Path(tempfile.mkdtemp(prefix=key+'.building-', dir=destination.parent))
            genes, models, regions = parse_annotation(gff)
            unavailable_products=[]
            aliases = defaultdict(set); sequences = defaultdict(list)
            for mid, item in models.items():
                for alias in item['aliases']:
                    if alias:
                        aliases[alias].add(mid)
            supplied = {};product_origins = {}
            if row.get('cds'):
                for header, coding in fasta_records(row['cds']):
                    code = re.search(r'\[transl_table=(\d+)\]', header)
                    if code and code[1] != '1':
                        continue
                    if len(coding)%3 or not set(coding)<=set('ACGT'):
                        continue
                    seq = translate(coding).removesuffix('*')
                    if not seq or '*' in seq:
                        continue
                    match = re.search(r'\[protein_id=([^\]]+)\]', header)
                    identifiers = {match[1]} if match else set(ACCESSION.findall(header))
                    for pid in identifiers:
                        for mid in aliases.get(pid, ()):
                            if mid in supplied and supplied[mid] != seq:
                                raise ValueError('Conflicting CDS products for transcript: '+mid)
                            supplied[mid] = seq;product_origins[mid] = 'translated_supplied_cds'
            if row.get('protein'):
                seen = set()
                for header, seq in fasta_records(row['protein']):
                    pid = header.split()[0]
                    if pid in seen:
                        raise ValueError('Duplicate protein FASTA identifier: '+pid)
                    seen.add(pid)
                    if not seq or not set(seq) <= AA:
                        continue
                    for mid in aliases.get(pid, ()):
                        supplied[mid] = seq;product_origins[mid] = 'supplied_annotation_protein'
            genome = None
            # A supplied proteome is mapped by exact annotation accession. Its
            # unmapped transcripts remain explicit, rather than using names.
            if models and use_genome:
                genome = Genome(fasta, tmp/'genome', runner)
            try:
                if genome:
                    for gene in genes.values():
                        if gene['end'] > genome.length(gene['contig']):
                            raise ValueError('GFF3/FASTA contig or length mismatch')
                    regions = {c:(0,v[0]) for c,v in genome.index.items()}
                for mid, item in models.items():
                    model = item['model']
                    if row.get('protein') or row.get('cds'):
                        seq, status = supplied.get(mid, ''), product_origins.get(mid, 'no_supplied_product')
                    else:
                        _, seq, status, issues = coding_sequence(model, genome)
                    if not seq or '*' in seq or not set(seq) <= AA or status == 'disrupted':
                        unavailable_products.append(dict(model_id=mid,gene_id=model.gene,reason=status if status=='disrupted' else 'no_usable_sequence',aliases=item['aliases']))
                        continue
                    sequences[model.gene].append(dict(sequence=seq, model_id=mid, aliases=item['aliases'], status=status))
            finally:
                if genome:
                    genome.close()
            records, sequence_map, alias_map, geneids = [], defaultdict(list), defaultdict(list), defaultdict(list)
            all_sequences={};product_count=0
            product_fields=['search_id','sequence_sha256','protein_ids','model_id','transcript_ids','gene_id','ncbi_gene_id','gene_name','contig','gene_start','gene_end','start','end','strand','protein_aa','cds','stop_codons','origin','exception','aliases']
            jsonl=(tmp/'products.jsonl').open('w');tsv=(tmp/'products.tsv').open('w',newline='')
            product_writer=csv.DictWriter(tsv,product_fields,delimiter='\t');product_writer.writeheader()
            for gid, gene in genes.items():
                gene['protein_id'] = ''; gene['protein_length'] = 0
                for product in sequences[gid]:
                    digest = hashlib.sha256(product['sequence'].encode()).hexdigest()
                    sid='A'+digest[:24];all_sequences[sid]=product['sequence']
                    item=models[product['model_id']];model=item['model']
                    transcript_ids=sorted({a for a in product['aliases'] if re.fullmatch(r'[NXY][MR]_\d+\.\d+',a)})
                    record=dict(search_id=sid,sequence_sha256=digest,protein_ids=item['protein_ids'],model_id=product['model_id'],
                        transcript_ids=transcript_ids,gene_id=gid,ncbi_gene_id=gene['ncbi_gene_id'],gene_name=gene['name'],
                        contig=gene['contig'],gene_start=gene['start'],gene_end=gene['end'],start=model.start,end=model.end,strand=model.strand,
                        protein_aa=len(product['sequence']),cds=sorted(model.cds),stop_codons=sorted(model.stops),origin=product['status'],
                        exception=model.attrs.get('exception',''),aliases=product['aliases'])
                    jsonl.write(json.dumps(record,separators=(',',':'))+'\n')
                    product_writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in record.items()})
                    product_count+=1
                    sequence_map[digest].append(gid)
                    for alias in product['aliases']:
                        if alias:
                            alias_map[alias].append([gid, digest])
                if gene['ncbi_gene_id']:
                    geneids[gene['ncbi_gene_id']].append(gid)
                if sequences[gid]:
                    best = min(sequences[gid], key=lambda p:(-len(p['sequence']), p['model_id']))
                    pid = 'G'+token(gid)[:24]
                    gene.update(protein_id=pid, protein_length=len(best['sequence']), protein_origin=best['status'])
                    records.append((pid, best['sequence']))
            jsonl.close();tsv.close()
            write_fasta(tmp/'proteins.faa', sorted(records))
            write_fasta(tmp/'all_products.faa',sorted(all_sequences.items()))
            save(tmp/'unavailable_products.json',unavailable_products)
            payload = dict(genes=list(genes.values()), aliases=dict(alias_map),
                sequence_genes={k:sorted(set(v)) for k,v in sequence_map.items()}, ncbi_gene_ids=dict(geneids),
                contigs=regions, annotation_genes=len(genes), coding_genes=sum(g['coding'] for g in genes.values()),
                protein_genes=len(records),all_product_models=product_count,unique_product_sequences=len(all_sequences),
                unavailable_product_models=len(unavailable_products),inputs=inputs,
                fasta_read=use_genome,contig_lengths_available=len(regions))
            save(tmp/'catalog.json', payload)
            # Genome copies/indexes are scratch data, not a persistent proteome.
            if (tmp/'genome').exists():
                shutil.rmtree(tmp/'genome')
            for p in input_paths:
                if (p.stat().st_size,p.stat().st_mtime_ns)!=initial_stats[str(p)]:
                    raise ValueError('Source dataset changed during catalog creation: '+str(p))
            seal_cache(tmp, identity);tmp.rename(destination)
    for p in input_paths:
        if (p.stat().st_size,p.stat().st_mtime_ns)!=initial_stats[str(p)]:
            raise ValueError('Source dataset changed while checking cache: '+str(p))
    payload = json.loads((destination/'catalog.json').read_text())
    return dict(payload, cache=str(destination), key=key, genome_id=row['genome_id'],
                species=row.get('species',''), assembly=row.get('assembly_accession',''),
                inputs=inputs, cache_reused=reused, elapsed_seconds=time.perf_counter()-t0,
                all_products_fasta=str(destination/'all_products.faa'),products_tsv=str(destination/'products.tsv'),
                products_jsonl=str(destination/'products.jsonl'))


def iter_products(catalog):
    """Stream every usable annotated product and binding, without loading an MSA.

    Several rows may share search_id: identical proteins are searched once but
    every transcript/gene/coordinate association remains a separate record.
    """
    with Path(catalog['products_jsonl']).open() as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


class GeneIndex:
    def __init__(self, catalog):
        self.genes = {g['gene_id']:g for g in catalog['genes']}
        self.starts, self.ends = {}, {}
        groups = defaultdict(list)
        for gene in self.genes.values():
            if gene['coding']:
                groups[gene['contig']].append(gene)
        for contig, values in groups.items():
            by_start = sorted(values, key=lambda g:(g['start'],g['end'],g['gene_id']))
            by_end = sorted(values, key=lambda g:(g['end'],g['start'],g['gene_id']))
            self.starts[contig] = ([g['start'] for g in by_start],by_start)
            self.ends[contig] = ([g['end'] for g in by_end],by_end)

    def neighbors(self, contig, start, end, strand, n, cap, exclude=()):
        excluded = set(exclude); result = []; limits = {}
        # All overlapping genes are omitted from neighbour anchors, including
        # alternative target models and nested/antisense genes.
        for side in ['left','right']:
            positions, values = (self.ends if side=='left' else self.starts).get(contig,([],[]))
            i = bisect_right(positions,start)-1 if side=='left' else bisect_left(positions,end)
            count=0; limited='contig_or_annotation_end'
            while 0<=i<len(values) and count<n:
                gene=values[i];i += -1 if side=='left' else 1
                if gene['gene_id'] in excluded:
                    continue
                if (side=='left' and gene['start']<max(0,start-cap)) or (side=='right' and gene['end']>end+cap):
                    limited='distance_cap';break
                count+=1
                result.append(dict(gene,side=side,side_rank=count,
                    target_side=('upstream' if side=='left' else 'downstream') if strand=='+' else ('downstream' if side=='left' else 'upstream'),
                    relative_strand='+' if gene['strand']==strand else '-',
                    distance_bp=start-gene['end'] if side=='left' else gene['start']-end))
            limits[side]='enough_genes' if count==n else limited
        return dict(neighbors=result,limits=limits,
                    window_start=min([start]+[g['start'] for g in result]),
                    window_end=max([end]+[g['end'] for g in result]))


def origin_tokens(member):
    identifiers=set()
    for k in ['query','protein','transcript']:
        value=member.get(k,'')
        identifiers.update(ACCESSION.findall(value))
        if value:
            identifiers.add(value)
    return identifiers


def resolve_origins(queries, sources):
    records=[];mapping=[]
    for query in queries:
        digest=hashlib.sha256(query['sequence'].encode()).hexdigest()
        for member in query['members']:
            strong=[];identical=[];mismatches=[]
            for cat in sources:
                explicit=member.get('genome_id') or member.get('source_genome_id')
                accession=member.get('assembly_accession')
                if explicit and explicit!=cat['genome_id']:
                    continue
                if accession and accession not in cat['genome_id'] and accession!=cat['assembly']:
                    continue
                if member.get('species') and cat['species'] and member['species']!=cat['species']:
                    continue
                ids=origin_tokens(member)
                linked={gid for name in ids for gid,seqhash in cat['aliases'].get(name,[]) if seqhash==digest}
                stale={gid for name in ids for gid,seqhash in cat['aliases'].get(name,[]) if seqhash!=digest}
                # GeneID is an origin hint only when an actual product sequence
                # also matches. A name/GeneID alone never validates the query.
                same=set(cat['sequence_genes'].get(digest,[]))
                if member.get('gene_id'):
                    hinted=set(cat['ncbi_gene_ids'].get(str(member['gene_id']),[]))
                    linked |= hinted & same
                strong.extend((cat['genome_id'],gid) for gid in linked)
                identical.extend((cat['genome_id'],gid) for gid in same)
                mismatches.extend((cat['genome_id'],gid) for gid in stale-linked)
            selected=sorted(set(strong or ([] if mismatches else identical)))
            status=('identifier_and_sequence' if strong and len(selected)==1 else
                    'multiple_source_loci' if len(selected)>1 else
                    'exact_sequence_only' if selected else
                    'identifier_sequence_mismatch' if mismatches else 'source_locus_unresolved')
            mapping.append(dict(query_id=query['id'],query=member['query'],status=status,
                                sources=[dict(genome_id=g,gene_id=gid) for g,gid in selected]))
            records.extend(dict(genome_id=g,gene_id=gid,query_id=query['id'],query=member['query'],resolution=status) for g,gid in selected)
    grouped={}
    for r in records:
        k=(r['genome_id'],r['gene_id'])
        grouped.setdefault(k,dict(genome_id=k[0],gene_id=k[1],queries=[],resolutions=[]))
        grouped[k]['queries'].append(r['query']);grouped[k]['resolutions'].append(r['resolution'])
    for r in grouped.values():
        r['queries']=sorted(set(r['queries']));r['resolutions']=sorted(set(r['resolutions']))
        r['reference_id']='R'+token([r['genome_id'],r['gene_id']])[:20]
    return list(grouped.values()), mapping


def summarize_hsps(hsps):
    qa, sa, chosen = [], [], []
    for h in sorted(hsps, key=lambda h:-h['bitscore']):
        q,s=(h['qstart'],h['qend']),(h['sstart'],h['send'])
        if any(max(q[0],a)<min(q[1],b) for a,b in qa) or any(max(s[0],a)<min(s[1],b) for a,b in sa):
            continue
        qa.append(q);sa.append(s);chosen.append(h)
    aligned=sum(h['length'] for h in chosen)
    return dict(bitscore=sum(h['bitscore'] for h in chosen),evalue=min(h['evalue'] for h in chosen),
        identity=sum(h['identity']*h['length'] for h in chosen)/aligned,
        coverage_query=sum(b-a for a,b in qa)/chosen[0]['qlen'],
        coverage_subject=sum(b-a for a,b in sa)/chosen[0]['slen'])


def match_panel(panel, target, cache_base, runner, versions, config, threads):
    """Search a small neighbour panel against the COMPLETE available proteome."""
    panel_identity=token(sorted(panel.items()))
    identity=dict(schema=SCHEMA,implementation_sha256=sha256(Path(__file__)),panel=panel_identity,target=target['key'],
                  tools=versions,thresholds={k:config[k] for k in ['min_coverage','min_identity','competing_ratio','evalue']})
    destination=cache_base/'matches'/token(identity)
    with cache_lock(destination.with_suffix('.lock')):
        reused=verified_cache(destination,identity)
        if not reused:
            tmp=Path(tempfile.mkdtemp(prefix=destination.name+'.building-',dir=destination.parent))
            write_fasta(tmp/'panel.faa',sorted(panel.items()))
            # The database also has its own content-addressed cache. New panels
            # do not rebuild the complete target protein database.
            db_identity=dict(schema=SCHEMA,catalog=target['key'],makeblastdb=versions['makeblastdb'])
            dbdir=cache_base/'databases'/token(db_identity)
            with cache_lock(dbdir.with_suffix('.lock')):
                if not verified_cache(dbdir,db_identity):
                    dbtmp=Path(tempfile.mkdtemp(prefix=dbdir.name+'.building-',dir=dbdir.parent))
                    runner.run(['makeblastdb','-in',Path(target['cache'])/'proteins.faa','-dbtype','prot','-out',dbtmp/'db'],dbtmp/'makeblastdb.stdout')
                    seal_cache(dbtmp,db_identity);dbtmp.rename(dbdir)
            runner.run(['blastp','-query',tmp/'panel.faa','-db',dbdir/'db',
                '-num_threads',threads,'-evalue',config['evalue'],'-max_target_seqs',max(1,target['protein_genes']),
                '-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'],tmp/'blastp.tsv')
            groups=defaultdict(list)
            target_ids={g['protein_id'] for g in target['genes'] if g['protein_id']}
            for line in (tmp/'blastp.tsv').read_text().splitlines():
                f=line.split('\t')
                if len(f)!=12 or f[0] not in panel or f[1] not in target_ids:
                    raise ValueError('Unexpected neighbour BLASTP output')
                qs,qe,ss,se=map(int,[f[4],f[5],f[6],f[7]])
                if not 1<=qs<=qe<=len(panel[f[0]]) or not 1<=ss<=se<=int(f[11]) or int(f[10])!=len(panel[f[0]]):
                    raise ValueError('Neighbour BLASTP lengths/coordinates differ')
                groups[f[0],f[1]].append(dict(qstart=qs-1,qend=qe,sstart=ss-1,send=se,
                    length=int(f[3]),identity=float(f[2])/100,bitscore=float(f[9]),evalue=float(f[8]),qlen=int(f[10]),slen=int(f[11])))
            by_query=defaultdict(list)
            for (qid,sid),hsps in groups.items():
                hit=dict(query_protein=qid,target_protein=sid,**summarize_hsps(hsps))
                if min(hit['coverage_query'],hit['coverage_subject'])>=config['min_coverage'] and hit['identity']>=config['min_identity']:
                    by_query[qid].append(hit)
            result={}
            for qid in panel:
                hits=sorted(by_query[qid],key=lambda h:(-h['bitscore'],h['evalue'],h['target_protein']))
                near=[h for h in hits if h['bitscore']>=hits[0]['bitscore']*config['competing_ratio']] if hits else []
                result[qid]=dict(status='unique_best_match' if len(near)==1 else 'ambiguous_matches' if near else 'no_match',
                                 best_matches=near,passing_matches=len(hits))
            save(tmp/'matches.json',result);seal_cache(tmp,identity);tmp.rename(destination)
    return dict(matches=json.loads((destination/'matches.json').read_text()),cache=str(destination),cache_reused=reused,
                panel_proteins=len(panel),target_proteins=target['protein_genes'])


def compare_context(reference, source_neighbors, target_neighbors, mapping):
    by_target={g['protein_id']:g for g in target_neighbors if g['protein_id']}
    # Source protein IDs are panel IDs with a source-catalog prefix, so genes
    # with identical symbols in different assemblies cannot collide.
    proposed=[];ambiguous=[]
    for gene in source_neighbors:
        qid=gene.get('panel_protein_id','')
        item=mapping.get(qid,dict(status='no_protein',best_matches=[]))
        local=[h for h in item['best_matches'] if h['target_protein'] in by_target]
        if item['status']!='unique_best_match':
            if local:
                ambiguous.append(dict(reference_gene=gene['gene_id'],matches=local))
            continue
        if local:
            h=local[0];other=by_target[h['target_protein']]
            proposed.append(dict(reference_gene=gene['gene_id'],reference_name=gene['name'],
                target_gene=other['gene_id'],target_name=other['name'],reference_side=gene['target_side'],
                target_side=other['target_side'],reference_start=gene['start'],target_start=other['start'],
                reference_relative_strand=gene['relative_strand'],target_relative_strand=other['relative_strand'],**h))
    multiplicity=Counter(m['target_gene'] for m in proposed)
    anchors=[m for m in proposed if multiplicity[m['target_gene']]==1]
    collisions=[m for m in proposed if multiplicity[m['target_gene']]>1]
    def position(m,which):
        strand=reference['strand'] if which=='reference' else reference['_target_strand']
        return m[which+'_start']*(1 if strand=='+' else -1)
    pairs=[(a,b) for i,a in enumerate(anchors) for b in anchors[i+1:]]
    order=sum((position(a,'reference')<position(b,'reference'))==(position(a,'target')<position(b,'target')) for a,b in pairs)/len(pairs) if pairs else None
    orientation=sum(m['reference_relative_strand']==m['target_relative_strand'] for m in anchors)/len(anchors) if anchors else None
    left=sum(m['reference_side']==m['target_side']=='upstream' for m in anchors)
    right=sum(m['reference_side']==m['target_side']=='downstream' for m in anchors)
    status=('insufficient_anchors' if len(anchors)<2 else
            'order_or_orientation_differs' if order!=1 or orientation!=1 else
            'two_sided_context_support' if left and right else 'one_sided_context_support')
    return dict(reference_id=reference['reference_id'],status=status,shared_anchors=len(anchors),
                left_anchors=left,right_anchors=right,order_concordance=order,orientation_concordance=orientation,
                matches=anchors,ambiguous_matches=ambiguous,colliding_matches=collisions,
                flags=['tentative_protein_homology_not_orthology']+(['competing_global_matches'] if ambiguous else [])+
                      (['several_reference_genes_match_one_target'] if collisions else []))


def run_early_synteny(data, *, genomes_manifest, source_genomes_manifest=None,
                      cache_dir, out_dir, threads=4, config=None, runner=None):
    """Return one self-contained JSON payload; do not edit a stage1 checkpoint.

    Called automatically before write_outputs/seal. A saved checkpoint can be
    enriched into a NEW checkpoint without another target-discovery run.
    """
    started=time.perf_counter();out=Path(out_dir);out.mkdir(parents=True,exist_ok=False)
    cache=Path(cache_dir).resolve()
    if cache==out.resolve() or cache in out.resolve().parents or out.resolve() in cache.parents:
        raise ValueError('Cache and evidence output must be separate directories')
    cache.mkdir(parents=True,exist_ok=True)
    cfg=dict(neighbors=5,max_neighbors=20,max_distance=1000000,min_coverage=.5,min_identity=.3,competing_ratio=.95,evalue=1e-5)
    cfg.update(config or {})
    if not (0<cfg['neighbors']<=cfg['max_neighbors'] and cfg['max_distance']>0):
        raise ValueError('Require 0 < neighbors <= max_neighbors and positive max_distance')
    if not (0<cfg['min_coverage']<=1 and 0<cfg['min_identity']<=1 and 0<cfg['competing_ratio']<=1 and cfg['evalue']>0):
        raise ValueError('Invalid preliminary homology thresholds')
    runner=runner or Runner(out)
    initial_commands=len(runner.commands)
    result=dict(schema=SCHEMA,status='running',parameters=cfg,source_resolution=[],datasets=[],references=[],loci=[],errors=[],commands=[])
    target_rows=datasets(genomes_manifest)
    source_rows=datasets(source_genomes_manifest) if source_genomes_manifest else target_rows
    # One physical dataset can be both a source reference and a target. Read it
    # once, while retaining dataset IDs used by the caller's stage1 loci.
    all_rows={}
    for r in target_rows+source_rows:
        gid=r['genome_id']
        if gid in all_rows:
            if any(r[k] and all_rows[gid][k] and r[k]!=all_rows[gid][k] for k in ['fasta','gff3','protein','cds','annotation_provenance']):
                raise ValueError('Same genome_id refers to conflicting source/target datasets: '+gid)
            all_rows[gid].update({k:v for k,v in r.items() if v})
        else:
            all_rows[gid]=dict(r)
    source_ids={r['genome_id'] for r in source_rows};catalogs={};indexes={};dataset_states={}
    known={g['genome_id']:g for g in data['genomes']}
    if data['config'].get('query_type')=='dna':
        result['status']='coding_reference_unavailable'
    else:
        for gid,row in all_rows.items():
            receipt=dict(genome_id=gid,role='source_and_target' if gid in source_ids and gid in known else 'source' if gid in source_ids else 'target')
            try:
                if not row['gff3']:
                    receipt['status']='annotation_unavailable'
                else:
                    print('[early synteny] prepare '+gid,flush=True)
                    old=known.get(gid,{})
                    digest=None
                    if old.get('fasta_sha256') and not (row.get('protein') or row.get('cds')):
                        stat=Path(row['fasta']).stat()
                        same_unchanged=(row['fasta']==old.get('fasta') and
                            stat.st_size==old.get('fasta_bytes') and stat.st_mtime_ns==old.get('fasta_mtime_ns'))
                        digest=old['fasta_sha256'] if same_unchanged else sha256(row['fasta'])
                        if digest!=old['fasta_sha256']:
                            raise ValueError('Target FASTA differs from the saved discovery input')
                    cat=build_catalog(row,cache,runner,digest)
                    catalogs[gid]=cat;indexes[gid]=GeneIndex(cat)
                    receipt.update(status='ready' if cat['protein_genes'] else 'protein_catalog_empty',
                        catalog_key=cat['key'],cache_reused=cat['cache_reused'],elapsed_seconds=cat['elapsed_seconds'],
                        annotation_genes=cat['annotation_genes'],coding_genes=cat['coding_genes'],protein_genes=cat['protein_genes'],inputs=cat['inputs'],
                        protein_strategy='longest_usable_product_per_gene; all usable products retained for source resolution')
            except (OSError,ValueError,RuntimeError,EOFError) as error:
                receipt.update(status='dataset_error',error=str(error));result['errors'].append(dict(genome_id=gid,error=str(error),stage='catalog'))
            dataset_states[gid]=receipt['status'];result['datasets'].append(receipt)
        sources=[cat for gid,cat in catalogs.items() if gid in source_ids]
        origins,resolution=resolve_origins(data['queries'],sources)
        result['source_resolution']=resolution
        target_genes=defaultdict(set)
        for ref in origins:
            target_genes[ref['genome_id']].add(ref['gene_id'])
        panel={};reference_contexts={}
        for ref in origins:
            gid=ref['genome_id'];cat=catalogs[gid];gene=indexes[gid].genes[ref['gene_id']]
            ref.update({k:gene[k] for k in ['name','contig','start','end','strand']})
            context=indexes[gid].neighbors(gene['contig'],gene['start'],gene['end'],gene['strand'],cfg['max_neighbors'],cfg['max_distance'],target_genes[gid])
            # Protein sequences are loaded only for the small neighbourhood panel.
            required={g['protein_id'] for g in context['neighbors'] if g['protein_id']}
            protein_sequences={header.split()[0]:seq for header,seq in fasta_records(Path(cat['cache'])/'proteins.faa') if header.split()[0] in required}
            for neighbor in context['neighbors']:
                pid=neighbor['protein_id']
                neighbor['panel_protein_id']='P'+hashlib.sha256(protein_sequences[pid].encode()).hexdigest()[:24] if pid else ''
                if pid:
                    panel[neighbor['panel_protein_id']]=protein_sequences[pid]
            ref.update(neighbors=context['neighbors'],limits=context['limits'],window_start=context['window_start'],window_end=context['window_end'])
            result['references'].append(ref);reference_contexts[ref['reference_id']]=context
        versions=tool_versions(['blastp','makeblastdb']);matches={}
        if panel and all(v['path'] for v in versions.values()):
            for row in target_rows:
                gid=row['genome_id']
                if gid not in catalogs or not catalogs[gid]['protein_genes']:
                    continue
                try:
                    print('[early synteny] match neighbour panel: '+gid,flush=True)
                    match=match_panel(panel,catalogs[gid],cache,runner,versions,cfg,threads)
                    matches[gid]=match['matches']
                    folder=out/'matches'/gid;folder.mkdir(parents=True)
                    for name in ['matches.json','blastp.tsv','blastp.tsv.log','cache.json']:
                        if (Path(match['cache'])/name).exists():shutil.copyfile(Path(match['cache'])/name,folder/name)
                    result['datasets'][[r['genome_id'] for r in result['datasets']].index(gid)].update(match_cache_reused=match['cache_reused'],panel_proteins=match['panel_proteins'])
                except (OSError,ValueError,RuntimeError,EOFError) as error:
                    result['errors'].append(dict(genome_id=gid,error=str(error),stage='protein_matching'))
        elif panel:
            result['errors'].append(dict(stage='protein_matching',error='Missing BLAST+ blastp/makeblastdb'))
        for locus in data['loci']:
            gid=locus['genome_id'];entry=dict(genome_id=gid,locus_id=locus['locus_id'],status='pending',
                neighbors_requested=cfg['neighbors'],neighbors_used=0,window_start=locus['start'],window_end=locus['end'],neighbors=[],comparisons=[],flags=[])
            if gid not in catalogs:
                entry['status']=dataset_states.get(gid,'dataset_unavailable')
            elif locus['contig'] not in catalogs[gid]['contigs']:
                entry['status']='annotation_contig_missing'
            elif catalogs[gid]['contigs'][locus['contig']][1]!=locus.get('contig_length',catalogs[gid]['contigs'][locus['contig']][1]):
                entry['status']='annotation_contig_length_mismatch'
            else:
                n=cfg['neighbors']
                while True:
                    context=indexes[gid].neighbors(locus['contig'],locus['start'],locus['end'],locus['strand'],n,cfg['max_distance'],target_genes[gid])
                    comparisons=[]
                    for ref in result['references']:
                        ref_neighbors=[g for g in ref['neighbors'] if g['side_rank']<=n]
                        comparisons.append(compare_context(dict(ref,_target_strand=locus['strand']),ref_neighbors,context['neighbors'],matches.get(gid,{})))
                    # Adaptive expansion is an annotation lookup, not another
                    # BLAST run. All reserve neighbours were searched in one panel.
                    if not comparisons or any(c['status']=='two_sided_context_support' for c in comparisons) or n>=cfg['max_neighbors']:
                        break
                    n=min(cfg['max_neighbors'],n*2)
                entry.update(neighbors_used=n,window_start=context['window_start'],window_end=context['window_end'],neighbors=context['neighbors'],comparisons=comparisons,limits=context['limits'])
                if not result['references']:
                    entry['status']='source_locus_unresolved'
                elif gid not in matches:
                    entry['status']='protein_matching_unavailable'
                elif any(c['status']=='two_sided_context_support' for c in comparisons):
                    entry['status']='two_sided_context_support'
                elif any(c['status']=='one_sided_context_support' for c in comparisons):
                    entry['status']='one_sided_context_support'
                elif any(c['shared_anchors']>=2 for c in comparisons):
                    entry['status']='order_or_orientation_differs'
                else:
                    entry['status']='insufficient_anchors'
                if n>cfg['neighbors']:entry['flags'].append('neighborhood_expanded_automatically')
                if any(not g['protein_id'] for g in context['neighbors']):entry['flags'].append('some_neighbors_lack_usable_proteins')
                if any(v!='enough_genes' for v in context['limits'].values()):entry['flags'].append('neighborhood_incomplete')
                entry['flags'].append('preliminary_context_evidence_not_orthology')
            result['loci'].append(entry)
    if not result['loci']:
        result['loci']=[dict(genome_id=l['genome_id'],locus_id=l['locus_id'],status=result['status'],neighbors_requested=cfg['neighbors'],neighbors_used=0,window_start=l['start'],window_end=l['end'],neighbors=[],comparisons=[],flags=[]) for l in data['loci']]
    result.update(status='completed_with_errors' if result['errors'] else 'completed',commands=runner.commands[initial_commands:],
                  elapsed_seconds=time.perf_counter()-started,tools=tool_versions(['blastp','makeblastdb']))
    save(out/'evidence.json',result)
    rows=[dict(genome_id=l['genome_id'],locus_id=l['locus_id'],status=l['status'],neighbors_requested=l['neighbors_requested'],neighbors_used=l['neighbors_used'],window_start=l['window_start'],window_end=l['window_end'],best_shared_anchors=max([c['shared_anchors'] for c in l['comparisons']] or [0]),flags=';'.join(l['flags'])) for l in result['loci']]
    write_tsv(out/'summary.tsv',rows,['genome_id','locus_id','status','neighbors_requested','neighbors_used','window_start','window_end','best_shared_anchors','flags'])
    write_fasta(out/'reference_neighbors.faa',sorted(panel.items()) if data['config'].get('query_type')!='dna' else [])
    return result


def main(argv=None):
    """Standalone enrichment of a sealed saved search; never mutates its files."""
    import argparse
    from .stage1 import verify
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',required=True,help='Existing sealed stage1 checkpoint')
    ap.add_argument('--genomes',required=True,help='Target datasets TSV: genome_id,fasta,gff3,optional protein/cds/species')
    ap.add_argument('--source-genomes',help='Optional source reference datasets; default: supplied target datasets')
    ap.add_argument('--cache',required=True)
    ap.add_argument('--out',required=True,help='New output directory, outside the saved checkpoint')
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--neighbors',type=int,default=5)
    ap.add_argument('--max-neighbors',type=int,default=20)
    ap.add_argument('--max-distance',type=int,default=1000000)
    args=ap.parse_args(argv)
    source,data=verify(args.results);out=Path(args.out).resolve()
    if source==out or source in out.parents or out in source.parents:
        raise ValueError('Enrichment output must be separate from source checkpoint')
    result=run_early_synteny(data,genomes_manifest=args.genomes,source_genomes_manifest=args.source_genomes,
        cache_dir=args.cache,out_dir=out,threads=args.threads,
        config=dict(neighbors=args.neighbors,max_neighbors=args.max_neighbors,max_distance=args.max_distance))
    print(json.dumps(dict(output=str(out),status=result['status'],reference_loci=len(result['references']),
        source_resolution=result['source_resolution'],locus_statuses=dict(Counter(l['status'] for l in result['loci'])),
        elapsed_seconds=result['elapsed_seconds'],errors=result['errors']),indent=2))
    return 1 if result['errors'] else 0


if __name__=='__main__':
    raise SystemExit(main())
