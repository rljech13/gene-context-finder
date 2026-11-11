"""Offline, locus-first review of existing BLAST evidence and early synteny."""
import json
from pathlib import Path

TEMPLATE = Path(__file__).with_name("stage1_review.txt").read_text().replace(
    "__CONTEXT_VIEW__", Path(__file__).with_name("stage1_context_view.txt").read_text())


def pack_report(data):
    """Intern repeated evidence objects without discarding a chain or isoform.

    The scientific checkpoint remains authoritative. This format is only an
    offline HTML transport: hydration restores the same report data, with shared
    read-only objects. Intern by full content, never by a possibly reused ID.
    """
    objects, seen = [], {}

    def visit(value, field=''):
        if isinstance(value, list):
            result = [visit(v, field + '[]') for v in value]
        elif isinstance(value, dict):
            result = {k: visit(v, k) for k, v in value.items()}
        else:
            return value
        if field in {'annotation_model', 'hsps[]', 'neighbors[]', 'comparisons[]',
                     'matches[]', 'colliding_matches[]', 'ambiguous_matches[]',
                     'display_only_genes[]', 'models', 'flags'} and value:
            key = json.dumps(result, sort_keys=True, separators=(',', ':'))
            if key not in seen:
                seen[key] = len(objects)
                objects.append(result)
            return {'$gcf_ref': seen[key]}
        return result

    # Receipts, commands and tool details remain in manifest/evidence, rather
    # than another copy of them in every standalone review page.
    view = {k: data[k] for k in ['review_id', 'config', 'filters', 'queries',
                                'genomes', 'loci', 'early_synteny', 'review_ranking',
                                'report_assets', 'copy_identities', 'copy_names_by_gene'] if k in data}
    return dict(schema='gcf_report_interned_v1', data=visit(view), objects=objects)


def write_report(root, data):
    from .gene_identity import identity, export_identities
    identities = {l['locus_id']: identity(l) for l in data['loci']}
    export_identities(root, data['loci'], identities, 'locus_id')
    references = data.get('early_synteny', {}).get('references', [])
    identities.update({r['reference_id']: identity(r) for r in references})
    copy_names = {}
    for record in list(data['loci']) + references:
        value = identity(record)
        if value['status'] == 'curated_locus_label':
            key = [record['genome_id'], record.get('annotation_gene_id') or record.get('gene_id')]
            copy_names[json.dumps(key, separators=(',', ':'))] = value['display_name']
    from .locus_ranking import rank_loci
    from .sequence import write_tsv
    ranking = rank_loci(data)
    fields = ['genome_id','locus_id','rank','score','strongest_source_weight','mean_source_weight','represented_sources','supplied_sources','discovery_rank']
    rows = sorted(ranking['loci'].values(), key=lambda r:(r['genome_id'],r['rank']))
    write_tsv(root/'review_ranking.tsv', rows, fields)
    manifest_path = root/'manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        manifest['review_ranking'] = dict(artifact='review_ranking.tsv', definition=ranking['definition'])
        manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    dna_path = root/'loci.fna'
    dna_loci = []
    if dna_path.exists():
        with dna_path.open() as handle:
            dna_loci = [line[1:].split()[0] for line in handle if line.startswith('>')]
    payload = json.dumps(pack_report(dict(data, review_ranking=ranking, copy_identities=identities, copy_names_by_gene=copy_names,
                         report_assets=dict(dna_loci=dna_loci))),
                         ensure_ascii=True, separators=(',', ':')).replace('<', chr(92)+'u003c').replace('&', chr(92)+'u0026')
    (root/'report.html').write_text(TEMPLATE.replace('__DATA__', payload))


def refresh(args):
    """Make a presentation revision while retaining the exact evidence/decision IDs."""
    import shutil
    from . import __version__
    from .sequence import sha256
    from .stage1 import verify, new_output, dump, now
    source, data = verify(args.results)
    target = Path(args.out).resolve()
    if target == source or source in target.parents or target in source.parents:
        raise ValueError('Report revision must be separate from the original checkpoint')
    root = new_output(target)
    original = json.loads((source/'seal.json').read_text())
    # Copy only checkpoint artifacts, never whole-genome scratch databases/workers.
    for name in original:
        destination = root/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, destination)
    if (source/'loci.decisions.tsv').exists():
        shutil.copyfile(source/'loci.decisions.tsv', root/'loci.decisions.tsv')
    manifest = json.loads((root/'manifest.json').read_text())
    manifest['presentation'] = dict(version=__version__, created=now(), source_checkpoint=str(source),
        source_evidence_sha256=sha256(source/'evidence.json'), retained_review_id=data['review_id'],
        renderer_sha256={name:sha256(Path(__file__).with_name(name)) for name in
                         ['stage1_report.py','stage1_review.txt','stage1_context_view.txt','locus_ranking.py','gene_identity.py']},
        note='Presentation revision only; all search evidence, filters, coordinates and review IDs retained.')
    dump(root/'manifest.json', manifest)
    write_report(root, data)
    dump(root/'seal.json', {name:sha256(root/name) for name in sorted(set(original) | {'review_ranking.tsv','gene_identities.tsv'})})
    print('Updated report: '+str(root/'report.html'))
    return 0
