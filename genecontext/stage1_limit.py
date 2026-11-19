"""Bound a one-gene review pool without treating omitted loci as rejected genes."""
from collections import Counter, defaultdict
from pathlib import Path
import shutil
from uuid import uuid4

MAX_LOCI_PER_ASSEMBLY = 20


def limit_candidates(data, maximum=MAX_LOCI_PER_ASSEMBLY):
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 1 <= maximum <= MAX_LOCI_PER_ASSEMBLY:
        raise ValueError('Candidate limit must be between 1 and 20 loci per gene per assembly')
    previous = data.get('candidate_limit')
    before = Counter(l['genome_id'] for l in data['loci'])
    if (previous and previous.get('scope') == 'one_gene_per_assembly'
            and previous['maximum_per_assembly'] == maximum and max(before.values(), default=0) <= maximum):
        return None
    from .locus_ranking import rank_loci
    ranking = rank_loci(data)['loci']
    pools = defaultdict(list)
    original = data['loci']
    for locus in original:
        pools[locus['genome_id']].append(locus)
    selected = [l for loci in pools.values()
                for l in sorted(loci, key=lambda l: ranking[l['locus_id']]['rank'])[:maximum]]
    kept = {l['locus_id'] for l in selected}
    data['loci'] = [l for l in original if l['locus_id'] in kept]
    after = Counter(l['genome_id'] for l in data['loci'])
    summary = dict(maximum_per_assembly=maximum, scope='one_gene_per_assembly',
                   policy='top_weighted_blast_review_rank_within_each_assembly',
                   source_pool_loci=len(original), discovered_loci=previous.get('discovered_loci',len(original)) if previous else len(original),
                   retained_loci=len(kept), omitted_from_pool=len(original)-len(kept),
                   by_genome={g:dict(source_pool=before[g],retained=after[g]) for g in sorted(before)},
                   audit='raw/candidate_limit.json', omission_is_biological_rejection=False)
    data['candidate_limit'] = summary
    data['config']['max_loci_per_assembly'] = maximum
    return dict(**summary, omitted_loci=[dict(locus_id=l['locus_id'],genome_id=l['genome_id'],
                review_rank=ranking[l['locus_id']]['rank']) for l in original if l['locus_id'] not in kept])


def apply_limit(args, root, data):
    if getattr(args, '_defer_candidate_limit', False):
        return
    maximum = getattr(args, 'max_loci_per_assembly', MAX_LOCI_PER_ASSEMBLY)
    audit = limit_candidates(data, maximum)
    from .stage1 import dump
    path = Path(root)/'raw/candidate_limit.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    if audit is not None or not path.exists():
        dump(path, audit if audit is not None else dict(data['candidate_limit'], inherited=True))
    if audit is not None and audit['omitted_from_pool']:
        subset_fasta(Path(root)/'loci.fna', {l['locus_id'] for l in data['loci']})


def subset_fasta(path, wanted, output=None):
    from .sequence import open_text
    if not path.exists():
        return
    target = output or path.with_suffix(path.suffix+'.limited')
    with open_text(path) as inp, target.open('w') as out:
        keep = False
        for line in inp:
            if line.startswith('>'):
                keep = line[1:].split()[0] in wanted
            if keep:
                out.write(line)
    if output is None:
        target.replace(path)


def inherit_parent_locator(source, root, data):
    """Retain an existing verified raw-evidence locator in derived reports."""
    if data.get('source_checkpoint'):
        parent = source/'_source'
        if not parent.is_dir():
            raise ValueError('Missing original checkpoint locator: '+str(parent))
        (root/'_source').symlink_to(parent.resolve(), target_is_directory=True)


def restrict_checkpoint(args):
    from .stage1 import verify, new_output, dump, write_outputs, now
    from .sequence import sha256
    source, data = verify(args.results)
    target = Path(args.out).resolve()
    if target == source or source in target.parents or target in source.parents:
        raise ValueError('Limited checkpoint must be separate from its source')
    root = new_output(target)
    (root/'raw').mkdir()
    old_review = data['review_id']
    source_receipt = dict(path=str(source),review_id=old_review,
                          evidence_sha256=sha256(source/'evidence.json'),seal_sha256=sha256(source/'seal.json'),
                          raw_base='_source/'+data.get('source_checkpoint',{}).get('raw_base','raw/'))
    data['source_checkpoint'] = source_receipt
    data['parent_review'] = old_review
    data['review_id'] = str(uuid4())
    data['created'] = now()
    apply_limit(args, root, data)
    # Keep one verified parent locator instead of copying its raw evidence tree.
    (root/'_source').symlink_to(source,target_is_directory=True)
    dump(root/'raw/source_checkpoint.json',source_receipt)
    shutil.copyfile(source/'queries.fasta',root/'queries.fasta')
    wanted={l['locus_id'] for l in data['loci']}
    subset_fasta(source/'loci.fna', wanted, root/'loci.fna')
    if data.get('early_synteny'):
        (root/'raw/early_synteny').mkdir()
        context=dict(data['early_synteny'],loci=[l['early_synteny'] for l in data['loci']])
        context['source_review_id']=old_review
        dump(root/'raw/early_synteny/evidence.json',context)
    write_outputs(root,data)
    print('STOP — open '+str(root/'report.html'))
    return 0
