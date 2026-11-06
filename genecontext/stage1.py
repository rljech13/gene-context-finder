"""One-gene, multi-reference BLAST discovery with a locus-first review stop.

Independent of the legacy model-first search. No predictor/reference downloads,
no orthology calls, and no automatic transition into context analysis.
"""
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
from uuid import uuid4

from . import __version__
from .comparison import table
from .locus_chaining import chain_hsps, group_loci, intersection_length, same_locus, summarise
from .pipeline import load_queries
from .stage1_proteins import discovery_rows, engine_for, search_annotated
from .sequence import read_gff, sha256, union_length, write_fasta, write_tsv
from .tools import Genome, Runner, tool_versions

OUTFMT = '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'
DECISIONS = ['review_id', 'genome_id', 'locus_id', 'decision', 'anchor_intervals',
             'window_start', 'window_end', 'notes', 'pair_decision', 'reason']
FILTERS = ['min_query_coverage', 'min_identity', 'max_evalue']
LOCUS_FIELDS = ['genome_id', 'rank_within_genome', 'locus_id', 'gene', 'contig', 'start', 'end', 'strand',
                'boundary_source', 'discovery_engine', 'annotation_gene_id', 'annotation_gene_name', 'query_variants', 'unique_sequences', 'support_species',
                'best_query', 'bitscore', 'coverage', 'identity', 'evalue', 'bitscore_per_query_unit',
                'pass_filters', 'filter_reasons', 'flags', 'window_start', 'window_end']


def now():
    return datetime.now(timezone.utc).isoformat()


def new_output(path):
    root = Path(path).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError('Output must be new or empty: '+str(root))
    root.mkdir(parents=True, exist_ok=True)
    return root


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def panel(args):
    originals = load_queries(args)
    metadata = {}
    if args.query_metadata:
        for row in table(Path(args.query_metadata)):
            key = row.get('query', '')
            if not key or key in metadata:
                raise ValueError('Query metadata requires unique query identifiers')
            if row.get('group') and row['group'] != args.gene:
                raise ValueError('This stage accepts one gene/search group; metadata group differs from --gene')
            metadata[key] = row
        if set(metadata) != {q['name'] for q in originals}:
            raise ValueError('Query metadata must describe every FASTA ID exactly once')
    unique, aliases = {}, {}
    for q in originals:
        # CDS route translates annotation-derived CDS before TBLASTn; DNA route
        # is strictly BLASTN, without automatic alphabet/type inference.
        sequence = q['protein'] if args.query_type in {'protein', 'cds'} else q['sequence']
        entry = dict(metadata.get(q['name'], {}), query=q['name'])
        if sequence not in unique:
            qid = f'q{len(unique)+1:04d}'
            unique[sequence] = dict(id=qid, sequence=sequence, members=[])
        value = unique[sequence]; value['members'].append(entry); aliases[q['name']] = value['id']
    return list(unique.values()), aliases


def parse_hits(path, queries, genome):
    result, seen = [], set()
    for number, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split('\t')
        if len(f) != 12 or f[0] not in queries or f[1] not in genome.index:
            raise ValueError(f'Malformed/unknown BLAST hit at line {number}')
        qstart, qend = int(f[4])-1, int(f[5])
        begin, end = int(f[6]), int(f[7])
        start, stop = min(begin, end)-1, max(begin, end)
        if not (0 <= qstart < qend <= len(queries[f[0]]['sequence'])) or not (0 <= start < stop <= genome.length(f[1])):
            raise ValueError('BLAST coordinates outside query/genome')
        if int(f[10]) != len(queries[f[0]]['sequence']) or int(f[11]) != genome.length(f[1]):
            raise ValueError('BLAST sequence lengths do not match inputs')
        if line in seen:
            continue
        seen.add(line)
        result.append(dict(hsp_id='h_'+hashlib.sha256(line.encode()).hexdigest()[:20],
            query=f[0], contig=f[1], start=start, end=stop, strand='+' if end>=begin else '-',
            qstart=qstart, qend=qend, alignment_length=int(f[3]), identity=float(f[2])/100,
            evalue=float(f[8]), bitscore=float(f[9]), qlen=int(f[10]), slen=int(f[11])))
    return result


def build_candidates(hsps, config, genome_id, queries):
    seeds, discarded = [], []
    for h in hsps:
        reasons = []
        if h['alignment_length'] < config['seed_min_length']: reasons.append('seed_too_short')
        if h['identity'] < config['seed_min_identity']: reasons.append('seed_low_identity')
        if reasons: discarded.append(dict(h, reasons=';'.join(reasons)))
        else: seeds.append(h)
    chains = chain_hsps(seeds, max_intron=config['max_intron'], max_overlap=config['chain_max_overlap'],
                        max_locus_span=config['max_locus_span'])
    loci = []
    for cluster in group_loci(chains, config['locus_shared_fraction']):
        first = cluster[0]
        start, end = min(c['start'] for c in cluster), max(c['end'] for c in cluster)
        token = json.dumps([genome_id, config['gene'], first['contig'], first['strand'],
                            sorted(c['chain_id'] for c in cluster)])
        locus_id = 'locus_'+hashlib.sha256(token.encode()).hexdigest()[:16]
        own_queries = sorted({c['query'] for c in cluster})
        members = [m for qid in own_queries for m in queries[qid]['members']]
        species = sorted({m['species'] for m in members if m.get('species')})
        limit = first['hsps'][0]['slen']
        loci.append(dict(genome_id=genome_id, locus_id=locus_id, gene=config['gene'],
            contig=first['contig'], start=start, end=end, strand=first['strand'],
            contig_length=limit, boundary_source='blast_chain_envelope',
            query_variants=len(members), unique_sequences=len(own_queries), support_species=len(species),
            species=species, query_ids=own_queries, chains=cluster, models=[],
            window_start=max(0,start-config['window']), window_end=min(limit,end+config['window']),
            flags=['gene_boundaries_unconfirmed']))
    # An envelope can contain/intersect another independent chain. Keep both
    # proposals and flag a possible fusion/fragmentation rather than hiding it.
    for i, a in enumerate(loci):
        for b in loci[i+1:]:
            if (a['contig'],a['strand']) == (b['contig'],b['strand']) and max(a['start'],b['start']) < min(a['end'],b['end']):
                a['flags'].append('overlapping_alternative_locus');b['flags'].append('overlapping_alternative_locus')
    return loci, chains, discarded


def refine_loci(loci, genome, runner, work, queries, config):
    """Refine every locus; no top-N cap and no removal on a failed alignment.

    Models must share aligned bases with the original seed blocks. Different
    query variants stay as separate models inside the locus. Conflicting model
    envelopes do not silently redefine the anchor or merge neighbouring loci.
    """
    if config['query_type'] == 'dna' or config['refine'] == 'none':
        for locus in loci: locus['flags'].append('refinement_not_requested')
        return
    if not shutil.which('miniprot'):
        for locus in loci: locus['flags'].append('miniprot_unavailable')
        return
    for locus in loci:
        folder = work/locus['locus_id'];folder.mkdir(parents=True)
        lo, hi = max(0,locus['start']-config['refine_padding']), min(genome.length(locus['contig']),locus['end']+config['refine_padding'])
        write_fasta(folder/'region.fna',[('region',genome.fetch(locus['contig'],lo,hi))])
        write_fasta(folder/'queries.faa',[(qid,queries[qid]['sequence']) for qid in locus['query_ids']])
        try:
            runner.run(['miniprot','-k','3','-M','0','-n','1','-l','3','-t',config['threads'],
                '-G',config['max_intron'],'-N',config['refine_max_models'],'--outn',config['refine_max_models'],'--outs','0.1',
                '--outc','0.1','--gff',folder/'region.fna',folder/'queries.faa'],folder/'models.gff3')
            models,_=read_gff(folder/'models.gff3',source='local_miniprot')
            retained=defaultdict(list)
            for m in models:
                target=m.attrs.get('Target','').split()
                if len(target)<3 or target[0] not in locus['query_ids'] or m.strand!=locus['strand']:
                    continue
                qid=target[0]; blocks=[(a+lo,b+lo) for a,b,_ in m.cds]
                seed=[block for c in locus['chains'] if c['query']==qid for block in c['blocks']]
                shared=intersection_length(blocks,seed)
                if not shared or shared/union_length(seed)<0.5:
                    continue
                retained[qid].append(dict(model_id=locus['locus_id']+'__'+m.id,
                    query=qid,start=m.start+lo,end=m.end+lo,
                    cds=[[a+lo,b+lo,p] for a,b,p in m.cds],
                    query_start=int(target[1])-1,query_end=int(target[2]),
                    coverage=(int(target[2])-int(target[1])+1)/len(queries[qid]['sequence']),
                    identity=float(m.attrs.get('Identity',0)),rank=int(m.attrs.get('Rank',1))))
            # Do not send hundreds of redundant low-ranking mapper alternatives
            # to the curator. Keep one best supported model per unique query;
            # the complete GFF3 (all alternatives) remains in raw evidence.
            locus['refinement_model_count']=sum(len(v) for v in retained.values())
            for values in retained.values():
                locus['models'].append(max(values,key=lambda m:(m['coverage']*m['identity'],-m['rank'])))
            locus['flags'].append('predicted_CDS_requires_review' if locus['models'] else 'CDS_not_recovered')
        except (OSError,ValueError,RuntimeError) as error:
            locus['flags'].append('refinement_error');locus['refinement_error']=str(error)


def consolidate_refined(loci, queries, config):
    """Use spliced aligned blocks to reconcile fragmentary BLAST proposals.

    Fixed model representatives are never enlarged by bridging chains. A chain
    that supports multiple anchors is retained as ambiguous evidence, without
    extending either anchor across two tandem copies. Unassigned chains remain
    independent BLAST-only candidates.
    """
    if not any(l['models'] for l in loci):return loci
    models={};chains={};origins=defaultdict(set)
    for l in loci:
        for c in l['chains']:
            chains[c['chain_id']]=c;origins[c['chain_id']].add(l['locus_id'])
        for m in l['models']:
            key=(m['query'],tuple(tuple(c) for c in m['cds']))
            models[key]=dict(m,contig=l['contig'],strand=l['strand'],
                blocks=[c[:2] for c in m['cds']],
                support_score=m['coverage']*m['identity']*len(queries[m['query']]['sequence']))
    clusters=[]
    for m in sorted(models.values(),key=lambda m:(-m['support_score'],m['model_id'])):
        matching=[g for g in clusters if same_locus(m,g[0],config['locus_shared_fraction'])]
        if matching:
            max(matching,key=lambda g:intersection_length(m['blocks'],g[0]['blocks'])).append(m)
        else:clusters.append([m])
    assigned=defaultdict(list);unassigned=[];ambiguous=set();split_anchors=set();weak_anchors=set()
    for c in chains.values():
        pieces=defaultdict(list);remaining=[];uncertain=set()
        for h in c['hsps']:
            support=[]
            for i,group in enumerate(clusters):
                representative=group[0]
                if (c['contig'],c['strand'])!=(representative['contig'],representative['strand']):continue
                shared=intersection_length([(h['start'],h['end'])],representative['blocks'])
                fraction=shared/(h['end']-h['start'])
                strong=fraction>=config['locus_shared_fraction']
                inside=representative['start']<=h['start']<h['end']<=representative['end']
                if strong or (inside and shared>0):support.append((shared,i,strong))
            if not support:remaining.append(h);continue
            shared,index,strong=max(support)
            pieces[index].append(h)
            if not strong:weak_anchors.add(index)
            if len(support)>1:
                uncertain.add(index);ambiguous.update(i for _,i,_ in support)
        if len(pieces)>1:split_anchors.update(pieces)
        for index,hsps in pieces.items():
            own=dict(c,hsps=hsps,blocks=[(h['start'],h['end']) for h in hsps],
                start=min(h['start'] for h in hsps),end=max(h['end'] for h in hsps),**summarise(hsps))
            if len(pieces)>1 or remaining:
                own['parent_chain_id']=c['chain_id']
                own['chain_id']=c['chain_id']+'__'+hashlib.sha256(';'.join(h['hsp_id'] for h in hsps).encode()).hexdigest()[:8]
                origins[own['chain_id']]=origins[c['chain_id']]
            if index in uncertain:own['ambiguous_support']=True
            assigned[index].append(own)
        if remaining:
            unassigned.append(dict(c,hsps=remaining,blocks=[(h['start'],h['end']) for h in remaining],
                start=min(h['start'] for h in remaining),end=max(h['end'] for h in remaining),**summarise(remaining)))
    result=[]
    for i,group in enumerate(clusters):
        if not assigned[i]:continue
        m=group[0];own=assigned[i]
        usable=[c for c in own if not c.get('ambiguous_support')] or own
        qids=sorted({c['query'] for c in usable})
        members=[v for qid in qids for v in queries[qid]['members']]
        source=next(l for l in loci if l['contig']==m['contig'] and l['strand']==m['strand'])
        # Boundaries come from an explicit representative model. Alternative
        # isoform bounds remain available for manual selection in the report.
        token=json.dumps([source['genome_id'],config['gene'],m['contig'],m['strand'],m['start'],m['end']])
        lid='locus_'+hashlib.sha256(token.encode()).hexdigest()[:16]
        species=sorted({v['species'] for v in members if v.get('species')})
        flags=['predicted_anchor_requires_review']
        if i in ambiguous:flags.append('ambiguous_chain_support')
        if i in split_anchors:flags.append('seed_chain_split_between_anchors')
        if i in weak_anchors:flags.append('partial_exon_overlap_support')
        overhang=max(30,round((m['end']-m['start'])*.05))
        if any(c['start']<m['start']-overhang or c['end']>m['end']+overhang for c in own):flags.append('BLAST_support_extends_beyond_anchor')
        # Models from repeated refinement windows can be identical. Collapse
        # them again and retain the best model for each query in this locus.
        by_query=defaultdict(list)
        for value in group:by_query[value['query']].append(value)
        kept=[max(v,key=lambda model:(model['coverage']*model['identity'],-model['rank'])) for v in by_query.values()]
        result.append(dict(genome_id=source['genome_id'],locus_id=lid,gene=config['gene'],
            contig=m['contig'],start=m['start'],end=m['end'],strand=m['strand'],
            contig_length=source['contig_length'],boundary_source='local_miniprot_representative',
            query_variants=len(members),unique_sequences=len(qids),support_species=len(species),species=species,
            query_ids=qids,chains=own,models=kept,flags=flags,
            merged_seed_loci=sorted({oid for c in own for oid in origins[c['chain_id']]}),
            window_start=max(0,m['start']-config['window']),window_end=min(source['contig_length'],m['end']+config['window'])))
    # Preserve unassigned evidence, without duplicating chains already attached
    # to model-backed loci or silently discarding a failed reconstruction.
    for group in group_loci(unassigned,config['locus_shared_fraction']):
        first=group[0];source=next(l for l in loci if l['locus_id'] in origins[first['chain_id']])
        start,end=min(c['start'] for c in group),max(c['end'] for c in group)
        members=[v for qid in {c['query'] for c in group} for v in queries[qid]['members']]
        qids=sorted({c['query'] for c in group});species=sorted({v['species'] for v in members if v.get('species')})
        lid='locus_'+hashlib.sha256(json.dumps([source['genome_id'],sorted(c['chain_id'] for c in group)]).encode()).hexdigest()[:16]
        result.append(dict(source,locus_id=lid,start=start,end=end,chains=group,models=[],query_ids=qids,
            boundary_source='blast_chain_envelope',query_variants=len(members),unique_sequences=len(qids),
            support_species=len(species),species=species,flags=['CDS_not_recovered','gene_boundaries_unconfirmed']+
                (['refinement_error'] if source.get('refinement_error') else []),
            window_start=max(0,start-config['window']),window_end=min(source['contig_length'],end+config['window'])))
    return sorted(result,key=lambda l:(l['contig'],l['start'],l['locus_id']))


def apply_filters(loci, filters, queries):
    for locus in loci:
        # Coverage and identity must pass in the SAME chain/query. Best values
        # from unrelated reference variants must not be combined to pass a gate.
        passing=[]
        for chain in locus['chains']:
            if chain.get('ambiguous_support'):continue
            if (chain['coverage']>=filters['min_query_coverage'] and chain['identity']>=filters['min_identity']
                    and chain['evalue']<=filters['max_evalue']):
                passing.append(chain)
        # Rank the strongest single-reference chain. Do not reward a locus for
        # many redundant isoforms, and do not change its score when filters move.
        pool=[c for c in locus['chains'] if not c.get('ambiguous_support')] or locus['chains']
        best=max(pool,key=lambda c:(c['bitscore'],c['coverage'],c['identity'],c['chain_id']))
        locus.update(best_query=best['query'],coverage=best['coverage'],identity=best['identity'],
            evalue=best['evalue'],bitscore=best['bitscore'],bitscore_per_query_unit=best['bitscore']/len(queries[best['query']]['sequence']),
            pass_filters=bool(passing),filter_reasons='' if passing else 'no_single_chain_passes_all_thresholds')
        locus['flags']=sorted(set(locus['flags']))
    groups=defaultdict(list)
    for locus in loci:groups[(locus['genome_id'],locus['gene'])].append(locus)
    order={key:i for i,key in enumerate(groups)}
    for group in groups.values():
        group.sort(key=lambda l:(-l['bitscore'],-l['coverage'],-l['identity'],l['contig'],l['start'],l['locus_id']))
        for rank,locus in enumerate(group,1):locus['rank_within_genome']=rank
    loci.sort(key=lambda l:(order[(l['genome_id'],l['gene'])],l['rank_within_genome']))


def summary_rows(data):
    rows=[]
    for genome in data['genomes']:
        loci=[l for l in data['loci'] if l['genome_id']==genome['genome_id']]
        rows.append(dict(genome_id=genome['genome_id'], status=genome['status'],
            discovery_engine=genome.get('discovery_engine','genomic'),
            searched_unique_products=genome.get('searched_unique_products',''),
            annotation_models=genome.get('all_product_models',''),unavailable_products=genome.get('unavailable_product_models',''),
            catalog_seconds=genome.get('catalog_seconds',''),blast_seconds=genome.get('blast_seconds',''),
            raw_hsps=genome.get('raw_hsps',0), loci=len(loci),passing_loci=sum(l['pass_filters'] for l in loci),
            refinement_errors=len(genome.get('refinement_errors',[])),
            control_expected=genome.get('control_expected',''),control_overlap=genome.get('control_overlap',''),
            error=genome.get('error','')))
    return rows


def write_outputs(root, data):
    from .stage1_report import write_report
    qmap={q['id']:q for q in data['queries']}
    apply_filters(data['loci'], data['filters'], qmap)
    dump(root/'evidence.json',data)
    rows=[]
    for locus in data['loci']:
        row={k:locus.get(k,'') for k in LOCUS_FIELDS};row['flags']=';'.join(locus['flags']);rows.append(row)
    write_tsv(root/'loci.tsv',rows,LOCUS_FIELDS)
    write_tsv(root/'passing_loci.tsv',[r for r in rows if r['pass_filters']],LOCUS_FIELDS)
    support=[]
    for locus in data['loci']:
        for c in locus['chains']:
            for member in qmap[c['query']]['members']:
                support.append(dict(locus_id=locus['locus_id'],genome_id=locus['genome_id'],
                    chain_id=c['chain_id'],query=member['query'],species=member.get('species',''),
                    gene_symbol=member.get('gene_symbol',''),transcript=member.get('transcript',''),
                    ambiguous_support=c.get('ambiguous_support',False),
                    query_length=len(qmap[c['query']]['sequence']),
                    query_aligned=union_length([(h['qstart'],h['qend']) for h in c['hsps']]),
                    query_intervals=';'.join(f"{h['qstart']}:{h['qend']}" for h in c['hsps']),
                    genomic_intervals=';'.join(f"{h['start']}:{h['end']}" for h in c['hsps'] if h.get('coordinate_space')!='protein'),
                    evidence_type=c.get('coordinate_space','genomic'),subject_protein=c.get('subject',''),
                    subject_length=c.get('subject_length',''),subject_coverage=c.get('subject_coverage',''),
                    subject_intervals=';'.join(f"{h['subject_start']}:{h['subject_end']}" for h in c['hsps'] if h.get('coordinate_space')=='protein'),
                    annotation_gene=locus.get('annotation_gene_id',''),annotation_model=c.get('annotation_model',{}).get('model_id',''),
                    contig=c['contig'],strand=c['strand'],
                    coverage=c['coverage'],identity=c['identity'],evalue=c['evalue'],bitscore=c['bitscore'],
                    hsp_ids=';'.join(h['hsp_id'] for h in c['hsps'])))
    write_tsv(root/'query_support.tsv',support,list(support[0]) if support else
              ['locus_id','genome_id','chain_id','query','species','coverage','identity','hsp_ids'])
    write_tsv(root/'genomes.summary.tsv',summary_rows(data),list(summary_rows(data)[0]))
    with open(root/'loci.bed','w') as handle:
        # Genome ID is a separate column: contig names need not be globally unique.
        for l in data['loci']:
            handle.write(f"{l['contig']}\t{l['start']}\t{l['end']}\t{l['locus_id']}\t0\t{l['strand']}\t{l['genome_id']}\n")
    decision_rows=[]
    for l in data['loci']:
        decision_rows.append(dict(review_id=data['review_id'],genome_id=l['genome_id'],
            locus_id=l['locus_id'],decision='',anchor_intervals=f"{l['start']}:{l['end']}",
            window_start='',window_end='',notes=''))
    for g in data['genomes']:
        if not any(l['genome_id']==g['genome_id'] for l in data['loci']):
            decision_rows.append(dict(review_id=data['review_id'],genome_id=g['genome_id'],
                locus_id='',decision='',anchor_intervals='',window_start='',window_end='',notes=''))
    write_tsv(root/'loci.decisions.tsv',decision_rows,DECISIONS)
    manifest=dict(version=__version__,stage='blast_locus_review',review_id=data['review_id'],
        created=now(),status='awaiting_review',config=data['config'],filters=data['filters'],
        query_provenance=data['query_provenance'],genomes=data['genomes'],tools=data['tools'],
        commands=data.get('commands',[]),parent_review=data.get('parent_review'),
        early_synteny=dict(status=data['early_synteny']['status'],artifact=data['early_synteny']['artifact']) if data.get('early_synteny') else None)
    dump(root/'manifest.json',manifest)
    write_report(root,data)
    files=[p for name in ['raw'] for p in (root/name).rglob('*') if p.is_file()]
    files += [root/name for name in ['evidence.json','manifest.json','queries.fasta','loci.tsv',
        'passing_loci.tsv','query_support.tsv','genomes.summary.tsv','loci.bed','report.html','loci.fna'] if (root/name).exists()]
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in sorted(files)})


def verify(root):
    root=Path(root).resolve()
    for name,digest in json.loads((root/'seal.json').read_text()).items():
        path=(root/name).resolve()
        if root not in path.parents or not path.is_file() or sha256(path)!=digest:
            raise ValueError('Stage 1 evidence changed: '+name)
    data=json.loads((root/'evidence.json').read_text())
    if data.get('schema')!='blast_loci_v1':raise ValueError('Not a locus-first stage 1 result')
    return root,data


def search(args):
    if not all(math.isfinite(v) for v in [args.max_evalue,args.search_evalue]) or args.max_evalue>args.search_evalue or args.max_evalue<=0 or args.search_evalue<=0:
        raise ValueError('Require 0 < --max-evalue <= --search-evalue')
    queries,aliases=panel(args);qmap={q['id']:q for q in queries}
    genomes=discovery_rows(args)
    engines={engine_for(row,args) for row in genomes}
    tools=tool_versions(['samtools','makeblastdb','blastp','tblastn','blastn','miniprot'])
    required=['makeblastdb',*engines]
    if engines-{'blastp'} or any(not (r.get('protein') or r.get('cds')) for r in genomes):required.append('samtools')
    if args.refine=='miniprot' and 'tblastn' in engines:required.append('miniprot')
    missing=[name for name in required if not tools[name]['path']]
    if missing:raise ValueError('Missing tools: '+', '.join(missing))
    root=new_output(args.out);(root/'raw').mkdir();runner=Runner(root)
    config={k:getattr(args,k) for k in ['gene','query_type','threads','search_evalue','seed_min_length',
        'seed_min_identity','max_intron','chain_max_overlap','max_locus_span','locus_shared_fraction',
        'window','refine','refine_padding','refine_max_models']}
    config['discovery']=getattr(args,'discovery','auto')
    data=dict(schema='blast_loci_v1',review_id=str(uuid4()),created=now(),queries=queries,
        query_provenance=dict(fasta=str(Path(args.query).resolve()),sha256=sha256(args.query),
            metadata=str(Path(args.query_metadata).resolve()) if args.query_metadata else None,
            metadata_sha256=sha256(args.query_metadata) if args.query_metadata else None,aliases=aliases),
        config=config,filters={k:getattr(args,k) for k in FILTERS},genomes=[],loci=[],tools=tools)
    write_fasta(root/'queries.fasta',[(q['id'],q['sequence']) for q in queries])
    sequences=[]
    for row in genomes:
        gid=row['genome_id'];record=dict(row,status='running');data['genomes'].append(record)
        work=root/'work'/gid;work.mkdir(parents=True);genome=None
        print('[stage 1] '+gid,file=sys.stderr,flush=True)
        try:
            record['discovery_engine']=engine_for(row,args)
            if record['discovery_engine']=='blastp':
                loci=search_annotated(row,args,config,qmap,root,work,runner,record)
                data['loci'].extend(loci)
                record['status']='completed' if loci else 'no_candidate'
                record['rescue_available']='Repeat this genome with --discovery genome if annotation support is absent or weak'
                continue
            initial_stat=Path(row['fasta']).stat()
            record.update(fasta_sha256=sha256(row['fasta']),fasta_bytes=initial_stat.st_size,
                          fasta_mtime_ns=initial_stat.st_mtime_ns)
            genome=Genome(row['fasta'],work,runner)
            database=work/'genome_db'
            runner.run(['makeblastdb','-in',genome.path,'-dbtype','nucl','-out',database],work/'makeblastdb.stdout')
            command=['blastn' if args.query_type=='dna' else 'tblastn','-query',root/'queries.fasta',
                '-db',database,'-evalue',args.search_evalue,'-num_threads',args.threads,
                '-max_target_seqs',max(1,len(genome.index)),'-outfmt',OUTFMT]
            if args.query_type=='dna':command+=['-task','blastn','-word_size','7','-strand','plus']
            output=root/'raw'/(gid+'.blast.tsv');runner.run(command,output)
            hsps=parse_hits(output,qmap,genome);record['raw_hsps']=len(hsps)
            loci,chains,discarded=build_candidates(hsps,config,gid,qmap)
            dump(root/'raw'/(gid+'.hsps.json'),dict(hsps=hsps,discarded_seeds=discarded))
            refine_loci(loci,genome,runner,work/'refinement',qmap,config)
            record['refinement_errors']=[dict(locus_id=l['locus_id'],error=l['refinement_error'])
                for l in loci if l.get('refinement_error')]
            raw_models=root/'raw'/(gid+'.refinement');raw_models.mkdir()
            for path in (work/'refinement').glob('*/models.gff3'):
                shutil.copyfile(path,raw_models/(path.parent.name+'.gff3'))
            loci=consolidate_refined(loci,qmap,config)
            own_sequences=[(l['locus_id'],genome.fetch(l['contig'],l['start'],l['end'],l['strand'])) for l in loci]
            final_stat=Path(row['fasta']).stat()
            if (initial_stat.st_size,initial_stat.st_mtime_ns)!=(final_stat.st_size,final_stat.st_mtime_ns):
                raise ValueError('Source FASTA changed while searching; results require a fresh run')
            sequences.extend(own_sequences)
            data['loci'].extend(loci);record['status']=('completed_with_refinement_errors' if record['refinement_errors']
                else 'completed' if loci else 'no_candidate')
            if row.get('expected_start') and row.get('expected_end'):
                record['control_expected']=f"{row['expected_start']}:{row['expected_end']}"
                record['control_overlap']=any(max(l['start'],int(row['expected_start']))<min(l['end'],int(row['expected_end']))
                    and l['strand']==row.get('expected_strand',l['strand']) for l in loci)
        except (OSError,ValueError,RuntimeError) as error:
            record.update(status='error',error=str(error));print(f'[stage 1] {gid}: {error}',file=sys.stderr)
        finally:
            if genome is not None:genome.close()
            dump(root/'progress.json',dict(updated=now(),genomes=data['genomes']))
    write_fasta(root/'loci.fna',sequences)
    data['commands']=runner.commands
    from .stage1_context import attach
    context_error=attach(args,root,data,runner)
    write_outputs(root,data)
    print('STOP — open '+str(root/'report.html'),file=sys.stderr)
    return 1 if context_error or any(g['status']=='error' or g.get('refinement_errors') for g in data['genomes']) else 0


def refilter(args):
    source,data=verify(args.results)
    for name in FILTERS:
        if getattr(args,name) is not None:data['filters'][name]=getattr(args,name)
    if not math.isfinite(data['filters']['max_evalue']) or not 0<data['filters']['max_evalue']<=data['config']['search_evalue']:
        raise ValueError('Review E-value must be positive and no looser than saved --search-evalue')
    root=new_output(args.out)
    if source in root.parents or root in source.parents:
        raise ValueError('Refilter output must be separate from the source')
    for name in ['queries.fasta','loci.fna']:
        if (source/name).exists():shutil.copyfile(source/name,root/name)
    shutil.copytree(source/'raw',root/'raw')
    data['parent_review']=data['review_id'];data['review_id']=str(uuid4())
    write_outputs(root,data)
    print('No BLAST rerun. New review: '+str(root/'report.html'),file=sys.stderr)
    return 0


def validate_decisions(data, rows):
    """Collect all decision errors before exporting any selected windows."""
    loci={l['locus_id']:l for l in data['loci']}
    genomes={g['genome_id']:g for g in data['genomes']}
    expected={(l['genome_id'],l['locus_id']) for l in data['loci']}
    expected|={(g['genome_id'],'') for g in data['genomes'] if not any(l['genome_id']==g['genome_id'] for l in data['loci'])}
    seen=set();approved=[];errors=[];by_genome=defaultdict(list)
    for number,row in enumerate(rows,2):
        key=(row.get('genome_id',''),row.get('locus_id',''))
        prefix=f'Row {number} ({key[0]} / {key[1] or "no candidate"}): '
        issues=[]
        if row.get('review_id')!=data['review_id']:issues.append('wrong review checkpoint')
        if key not in expected:issues.append('unknown genome/locus')
        if key in seen:issues.append('duplicate genome/locus')
        if issues:
            errors.extend(prefix+message for message in issues);continue
        seen.add(key);by_genome[key[0]].append(row)
        action=(row.get('decision') or '').strip();notes=(row.get('notes') or '').strip()
        reason=(row.get('reason') or '').strip();explanation=reason or notes
        pair_action=(row.get('pair_decision') or '').strip()
        if pair_action and pair_action not in {'inspect','exclude','unresolved'}:
            errors.append(prefix+'invalid pair decision')
        if action not in {'inspect','exclude','unresolved'}:
            errors.append(prefix+'pending or invalid decision');continue
        if action in {'exclude','unresolved'}:
            if not explanation:errors.append(prefix+'a short reason or note is required')
            continue
        try:
            if not key[1]:raise ValueError('cannot inspect a genome without a discovered locus')
            if genomes[key[0]]['status']=='error':raise ValueError('search failed; resolve the failed search before selecting a locus')
            l=loci[key[1]]
            intervals=(row.get('anchor_intervals') or '').strip() or f"{l['start']}:{l['end']}"
            anchors=[]
            for value in intervals.split(';'):
                try:a,b=map(int,value.split(':'))
                except ValueError:raise ValueError('anchor intervals must be start:end separated by semicolons')
                if not (0<=a<b<=l['contig_length']) or max(a,l['start'])>=min(b,l['end']):
                    raise ValueError('anchor must overlap the candidate and stay on its contig')
                anchors.append((a,b))
            if len(set(anchors))!=len(anchors):raise ValueError('duplicate manual anchor intervals')
            changed=anchors!=[(l['start'],l['end'])]
            if (changed or not l['pass_filters']) and not explanation:
                raise ValueError('boundary changes or below-filter selections need a short reason')
            start=(row.get('window_start') or '').strip();end=(row.get('window_end') or '').strip()
            custom=bool(start or end)
            if custom and (len(anchors)!=1 or not start or not end):raise ValueError('custom window requires both bounds and one anchor')
            if custom and not explanation:raise ValueError('manual window changes need a short reason')
            own=[]
            for i,(a,b) in enumerate(anchors,1):
                ws=int(start) if custom else max(0,a-data['config']['window'])
                we=int(end) if custom else min(l['contig_length'],b+data['config']['window'])
                if not (0<=ws<=a<b<=we<=l['contig_length']):raise ValueError('window must contain the anchor and stay on the contig')
                own.append(dict(genome_id=l['genome_id'],locus_id=l['locus_id'],region_id=l['locus_id']+(f'__part{i}' if len(anchors)>1 else ''),
                    gene=l['gene'],contig=l['contig'],strand=l['strand'],anchor_start=a,anchor_end=b,window_start=ws,window_end=we,
                    boundary_source='manual' if changed else l['boundary_source'],reason=reason,notes=notes))
            approved.extend(own)
        except (ValueError,TypeError) as error:errors.append(prefix+str(error))
    for gid,lid in sorted(expected-seen):errors.append(f'Missing decision: {gid} / {lid or "no candidate"}')
    pairs=[]
    for gid in genomes:
        own=by_genome[gid];declared={r.get('pair_decision','').strip() for r in own if r.get('pair_decision','').strip()}
        actions={r.get('decision','').strip() for r in own};selected=sum(r.get('decision')=='inspect' for r in own)
        if declared:
            if len(declared)!=1 or len(own)!=sum(bool(r.get('pair_decision','').strip()) for r in own):
                errors.append(gid+': inconsistent or incomplete pair decision metadata')
            else:
                decision=next(iter(declared))
                if decision=='inspect' and (not selected or not actions<={'inspect','exclude'}):errors.append(gid+': inspect pair must select loci and exclude the other loci')
                if decision in {'exclude','unresolved'} and actions!={decision}:errors.append(gid+': locus actions conflict with pair decision')
        decision=next(iter(declared)) if len(declared)==1 else ('inspect' if selected else next(iter(actions)) if len(actions)==1 else 'mixed_locus_decisions')
        pairs.append(dict(genome_id=gid,gene=data['config']['gene'],decision=decision,selected_loci=selected,
                          regions=sum(r['genome_id']==gid for r in approved)))
    return dict(valid=not errors,review_id=data['review_id'],errors=errors,pairs=pairs,regions=approved,
                selected_loci=sum(p['selected_loci'] for p in pairs),region_count=len(approved))


def approve(args):
    source,data=verify(args.results)
    result=validate_decisions(data,table(args.decisions))
    if getattr(args,'validate_only',False):
        print(json.dumps(result,indent=2))
        return 0 if result['valid'] else 2
    if not result['valid']:raise ValueError('Decision validation failed:\n'+'\n'.join(result['errors']))
    if not args.out:raise ValueError('--out is required unless --validate-only is used')
    root=new_output(args.out);approved=result['regions']
    fields=['genome_id','locus_id','region_id','gene','contig','strand','anchor_start','anchor_end','window_start','window_end','boundary_source','reason','notes']
    write_tsv(root/'approved_loci.tsv',approved,fields)
    with open(root/'synteny_windows.bed','w') as handle:
        for r in approved:handle.write(f"{r['contig']}\t{r['window_start']}\t{r['window_end']}\t{r['region_id']}\t0\t{r['strand']}\t{r['genome_id']}\n")
    shutil.copyfile(args.decisions,root/'decisions.snapshot.tsv');dump(root/'validation.json',result)
    dump(root/'manifest.json',dict(stage='approved_loci',review_id=data['review_id'],created=now(),source=str(source),
        source_evidence_sha256=sha256(source/'evidence.json'),decisions_sha256=sha256(args.decisions),regions=len(approved),
        note='Explicitly reviewed regions; not an orthology assertion. Context prediction has not run.'))
    return 0
