"""Annotation-first discovery; protein coordinates never masquerade as DNA hits."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import time

from .locus_chaining import chain_hsps
from .sequence import union_length, write_fasta


def discovery_rows(args):
    from .pipeline import load_genomes
    rows = load_genomes(args)
    from .early_synteny import datasets
    normalized={r['genome_id']:r for r in datasets(args.genomes)}
    for row in rows:row.update(normalized[row['genome_id']])
    if getattr(args, 'datasets', None):
        extra = {r['genome_id']: r for r in datasets(args.datasets)}
        for row in rows:
            for key, value in extra.get(row['genome_id'], {}).items():
                if value and row.get(key) and row[key] != value:
                    raise ValueError('Conflicting search/dataset field: '+row['genome_id']+' / '+key)
                if value:
                    row[key] = value
    return rows


def engine_for(row, args):
    if args.query_type == 'dna':
        if getattr(args, 'discovery', 'auto') == 'blastp':
            raise ValueError('BLASTP requires protein or CDS queries; DNA uses BLASTN')
        return 'blastn'
    mode = getattr(args, 'discovery', 'auto')
    return 'blastp' if mode == 'blastp' or (mode == 'auto' and row.get('gff3')) else 'tblastn'


def parse_protein_hits(path, queries, subjects):
    hits, seen = [], set()
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if line in seen:
            continue
        seen.add(line)
        f = line.split('\t')
        if len(f) != 12 or f[0] not in queries or f[1] not in subjects:
            raise ValueError(f'Unknown/malformed BLASTP hit, line {n}')
        qa, qb, sa, sb = int(f[4])-1, int(f[5]), int(f[6])-1, int(f[7])
        qlen, slen = len(queries[f[0]]['sequence']), subjects[f[1]]
        if not (0 <= qa < qb <= qlen and 0 <= sa < sb <= slen) or (int(f[10]), int(f[11])) != (qlen, slen):
            raise ValueError(f'BLASTP coordinates/lengths differ from input, line {n}')
        hits.append(dict(hsp_id='h_'+hashlib.sha256(line.encode()).hexdigest()[:20],
            query=f[0], contig=f[1], start=sa, end=sb, strand='+',
            qstart=qa, qend=qb, alignment_length=int(f[3]), identity=float(f[2])/100,
            evalue=float(f[8]), bitscore=float(f[9]), qlen=qlen, slen=slen))
    return hits


def search_annotated(row, args, config, queries, root, work, runner, record):
    """One annotated gene is one locus, including every matched transcript."""
    from .early_synteny import build_catalog, fasta_records, iter_products
    from .stage1 import OUTFMT, dump
    t0 = time.perf_counter()
    if not row.get('gff3'):
        raise ValueError('Annotation-first discovery requires GFF3 with gene/transcript/CDS relationships')
    cache = Path(getattr(args, 'neighbor_cache', None) or root.parent/'_genecontext_cache'/'early_synteny').resolve()
    catalog = build_catalog(row, cache, runner)
    fasta = Path(catalog['cache'])/'all_products.faa'
    # Public product export uses one sequence search ID with all gene/model
    # mappings. Do not use proteins.faa (one representative per neighbor gene).
    lengths = {pid: len(seq) for header, seq in fasta_records(fasta) for pid in [header.split()[0]]}
    if not lengths:
        raise ValueError('No usable annotation products; consider genomic rescue')
    genes = {g['gene_id']: g for g in catalog['genes']}
    record.update(discovery_engine='blastp', catalog_cache=catalog['cache'],
        catalog_key=catalog['key'], annotation_inputs=catalog['inputs'],
        catalog_reused=catalog['cache_reused'], catalog_seconds=catalog['elapsed_seconds'],
        searched_unique_products=len(lengths), annotation_genes=len(genes),
        all_product_models=catalog['all_product_models'], unavailable_product_models=catalog['unavailable_product_models'],
        fasta_read=catalog['fasta_read'],
        dna_extraction='deferred_until_locus_selection')
    database = work/'protein_db'
    runner.run(['makeblastdb', '-in', fasta, '-dbtype', 'prot', '-out', database], work/'makeblastdb.stdout')
    output = root/'raw'/(row['genome_id']+'.blast.tsv')
    blast_start = time.perf_counter()
    runner.run(['blastp', '-query', root/'queries.fasta', '-db', database,
        '-evalue', args.search_evalue, '-num_threads', args.threads,
        '-max_target_seqs', len(lengths), '-outfmt', OUTFMT], output)
    record['blast_seconds'] = time.perf_counter()-blast_start
    hits = parse_protein_hits(output, queries, lengths)
    hit_subjects={h['contig'] for h in hits}
    by_subject={pid:dict(search_id=pid,mappings=[]) for pid in hit_subjects}
    for product in iter_products(catalog):
        if product['search_id'] in hit_subjects:by_subject[product['search_id']]['mappings'].append(product)
    if any(not p['mappings'] for p in by_subject.values()):
        raise ValueError('BLAST subject has no annotation product mapping')
    seeds, discarded = [], []
    for h in hits:
        reasons = []
        if h['alignment_length'] < config['seed_min_length']: reasons.append('seed_too_short')
        if h['identity'] < config['seed_min_identity']: reasons.append('seed_low_identity')
        if reasons: discarded.append(dict(h, reasons=';'.join(reasons)))
        else: seeds.append(h)
    # Collinearity is evaluated on a single query/subject protein. Intron and
    # genomic-span limits are irrelevant here; subject length is a strict bound.
    longest = max(lengths.values())
    chains = chain_hsps(seeds, max_intron=longest,
        max_overlap=config['chain_max_overlap'], max_locus_span=longest)
    assigned = defaultdict(list)
    matched_products = set()
    for c in chains:
        pid = c['contig']; product = by_subject[pid]; matched_products.add(pid)
        for model in product['mappings']:
            gid = model['gene_id']
            if gid not in genes:
                raise ValueError('Product references unknown annotation gene: '+gid)
            gene = genes[gid]
            hsps = []
            for h in c['hsps']:
                item = {k:v for k,v in h.items() if k not in {'contig','strand','start','end'}}
                item.update(coordinate_space='protein', subject=pid,
                            subject_start=h['start'], subject_end=h['end'])
                hsps.append(item)
            own = dict(c, chain_id='chain_'+hashlib.sha256(json.dumps([c['chain_id'],gid,model['model_id']]).encode()).hexdigest()[:16],
                contig=gene['contig'], strand=gene['strand'], start=gene['start'], end=gene['end'],
                blocks=[b[:2] for b in model['cds']], hsps=hsps,
                coordinate_space='protein', subject=pid, subject_length=lengths[pid],
                subject_coverage=union_length([(h['start'],h['end']) for h in c['hsps']])/lengths[pid],
                annotation_model=model, shared_product_loci=sorted({m['gene_id'] for m in product['mappings']}))
            assigned[gid].append(own)
    loci = []
    for gid, own in assigned.items():
        gene = genes[gid]; qids = sorted({c['query'] for c in own})
        members = [m for qid in qids for m in queries[qid]['members']]
        species = sorted({m['species'] for m in members if m.get('species')})
        region = catalog['contigs'].get(gene['contig'])
        if not region or region[1] < gene['end']:
            raise ValueError('Missing/invalid contig extent for annotated locus: '+gene['contig'])
        lid = 'locus_'+hashlib.sha256(json.dumps([row['genome_id'], config['gene'], gid]).encode()).hexdigest()[:16]
        flags = ['annotation_gene_boundaries', 'annotation_completeness_not_guaranteed']
        if any(len(c['shared_product_loci'])>1 for c in own): flags.append('identical_product_at_multiple_gene_loci')
        loci.append(dict(genome_id=row['genome_id'], locus_id=lid, gene=config['gene'],
            annotation_gene_id=gid, annotation_gene_name=gene.get('name',gid),
            contig=gene['contig'], start=gene['start'], end=gene['end'], strand=gene['strand'],
            contig_length=region[1], boundary_source='annotated_gene', discovery_engine='blastp',
            query_variants=len(members), unique_sequences=len(qids), support_species=len(species),
            species=species, query_ids=qids, chains=own, models=[],
            window_start=max(0,gene['start']-config['window']),
            window_end=min(region[1],gene['end']+config['window']), flags=flags))
    folder = root/'raw'/(row['genome_id']+'.annotation'); folder.mkdir()
    dump(folder/'products.json', dict(products=[by_subject[p] for p in sorted(matched_products)]))
    dump(folder/'catalog_receipt.json', {k:catalog.get(k) for k in
        ['key','inputs','annotation_genes','coding_genes','protein_genes','all_product_models',
         'unique_product_sequences','unavailable_product_models','fasta_read','cache_reused','elapsed_seconds']})
    # Preserve explicit gaps in the product inventory for review/rescue.
    import shutil
    shutil.copyfile(Path(catalog['cache'])/'unavailable_products.json',folder/'unavailable_products.json')
    write_fasta(folder/'matched_products.faa', [(header,seq) for header,seq in fasta_records(fasta) if header.split()[0] in matched_products])
    dump(root/'raw'/(row['genome_id']+'.hsps.json'), dict(coordinate_space='protein',hsps=hits,discarded_seeds=discarded))
    record.update(raw_hsps=len(hits), refinement_errors=[], matched_annotated_genes=len(loci),
                  discovery_seconds=time.perf_counter()-t0)
    return loci
