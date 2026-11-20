"""Alignment-first STOP 2: candidate membership is distinct from approval.

The adapter reads sealed STOP 1 products or reuses a sealed legacy STOP 2 pool.
It never searches, runs Tiberius, modifies old checkpoints or transfers approval.
"""
from html import escape
import hashlib
import json
from pathlib import Path
from uuid import uuid4

from .alignment_blossom import BLOSSOM
from .sequence import read_fasta, sha256
from .stage1 import dump, new_output, now, verify

SCHEMA = 'gcf_alignment_first_v1'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def legacy_pool(source):
    from .stage2_choices import verify_checkpoint, constraints
    from .interactive_review import stage2_dataset
    source, data = verify_checkpoint(source)
    pool = stage2_dataset(source, data)
    rules = constraints(data)
    target_physical = {rules[m['model_id']]['physical_id'] for m in data['models']
                       if m['overlaps_anchor'] and m['strand'] == next(r['strand'] for r in data['regions'] if r['region_id'] == m['region_id'])}
    pool['sequences'] = [s for s in pool['sequences'] if s['reference'] or
                         rules[s['id']]['physical_id'] in target_physical]
    for s in pool['sequences']:
        if not s['reference']:
            s['metadata']['model_constraints'] = rules[s['id']]
    pool['selected'] = [sid for sid in pool['visible'] if any(s['id'] == sid for s in pool['sequences'])]
    return pool, dict(kind='legacy_stop2_reuse', checkpoint=str(source),
        source_review_id=data['review_id'], evidence_sha256=sha256(source/'evidence.json'),
        seal_sha256=sha256(source/'seal.json'), tiberius='available_before_model_choice_context_view_after_tree',
        copied_approvals=False)


def stage1_pool(source, args):
    """Use exact sealed BLASTP products; genomic rescue gaps remain explicit."""
    from .stage2 import selected_regions
    source, data = verify(source)
    regions, selection, _ = selected_regions(data, args)
    loci = {l['locus_id']: l for l in data['loci']}
    records, missing = [], []
    for region in regions:
        locus = loci[region['locus_id']]
        folder = source/'raw'/(region['genome_id']+'.annotation')
        if locus.get('discovery_engine') != 'blastp' or not (folder/'products.json').exists():
            missing.append(dict(locus_id=locus['locus_id'], reason='No sealed annotation protein products; requires targeted rescue after review.'))
            continue
        products = json.loads((folder/'products.json').read_text())['products']
        proteins = read_fasta(folder/'matched_products.faa')
        seen = set()
        for product in products:
            for m in product['mappings']:
                if m['gene_id'] != locus.get('annotation_gene_id') or m['model_id'] in seen:
                    continue
                seen.add(m['model_id'])
                seq = proteins[product['search_id']]
                if hashlib.sha256(seq.encode()).hexdigest() != m['sequence_sha256']:
                    raise ValueError('STOP 1 product identity differs: '+m['model_id'])
                mid = 'A'+digest([region['genome_id'], m['model_id'], m['sequence_sha256']])[:24]
                records.append(dict(id=mid, protein=seq, cds='', group=region['gene'], reference=False,
                    label=region['genome_id'].split('_GC')[0].replace('_', ' ')+' · '+m['model_id']+' · annotation',
                    metadata=dict(m, genome_id=region['genome_id'], locus_id=locus['locus_id'],
                                  stop1_comment=region.get('notes', ''), source='annotation',
                                  cds_status='not_exported_by_stop1', approval='not_reviewed')))
    for q in data['queries']:
        records.append(dict(id='reference_'+q['id'], protein=q['sequence'], cds='',
            group=data['config']['gene'], reference=True,
            label='Query · '+q['id']+' · '+', '.join(m.get('species','') for m in q['members']),
            metadata=dict(aliases=q['members'], cds_status='not_supplied')))
    return dict(sequences=records, selected=[s['id'] for s in records],
                visible=[s['id'] for s in records], colours=BLOSSOM), dict(
        kind='stop1_annotation_products', checkpoint=str(source), source_review_id=data['review_id'],
        evidence_sha256=sha256(source/'evidence.json'), seal_sha256=sha256(source/'seal.json'),
        selection=selection, missing_loci=missing, tiberius='not_run', copied_approvals=False)


def prepared_page(source, pool, title):
    """Re-render the sealed native payload in the current review order."""
    from .interactive_review import attach
    from .stage2_report import render_document

    seal = json.loads((source/'seal.json').read_text())
    report = source/'report.html'
    if seal.get('report.html') != sha256(report):
        raise ValueError('The prepared model report must match its evidence seal')
    text = report.read_text()
    offset = text.index('const DATA=')+len('const DATA=')
    payload, _ = json.JSONDecoder().raw_decode(text[offset:])
    if payload['review_id'] != pool['source_review_id']:
        raise ValueError('Prepared model report belongs to another review')
    if payload.get('annotation_context_external'):
        supplement = source/'annotation_context.json'
        if seal.get('annotation_context.json') != sha256(supplement):
            raise ValueError('Prepared annotation context must match its seal')
        payload['annotation_features'] = json.loads(supplement.read_text())['features']
        payload['annotation_context_external'] = False
    html = attach(render_document(payload), pool)
    html = html.replace('</head>', '<meta name="gcf-review-flow" content="sequential-v2"></head>', 1)
    scripts = []
    for name in ('alignment_review_steps.txt','synteny_review_page.txt','stop4_handoff.txt'):
        script = Path(__file__).with_name(name).read_text()
        scripts.append('<script>'+script+'</script>')
    from .synteny_reference import attach as attach_reference
    return attach_reference(html.replace('</body>', ''.join(scripts)+'</body>', 1))


def build(source, out, *, legacy=False, args=None, title='Candidate models'):
    pool, provenance = legacy_pool(source) if legacy else stage1_pool(source, args)
    if not any(not s['reference'] for s in pool['sequences']):
        raise ValueError('No exact candidate proteins available for alignment review')
    # A new review namespace isolates model membership from old approvals/drafts.
    pool.update(schema='gcf_interactive_review_v1', workflow=SCHEMA, stage=2, title=title,portable_available=False,
        review_id=str(uuid4()), source_review_id=provenance['source_review_id'],
        history=provenance, phylogeny=dict(method='ML, ModelFinder BIC, UFBoot 1000, seed 42',
            rooting='unrooted_without_validated_outgroup', interpretation='Candidate comparison, not orthology or mRNA approval'))
    pool['visible'] = list(dict.fromkeys(pool['visible']))
    pool['selected'] = list(dict.fromkeys(pool['selected']))
    root = new_output(out)
    dump(root/'review_data.json', pool)
    from .interactive_review import attach
    html = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+escape(title)+' · alignment-first STOP 2</title></head><body><main style="max-width:1400px;margin:24px auto;padding:16px;font-family:system-ui"><h1>'+escape(title)+' · STOP 2: alignment and gene tree</h1><p>STOP 1 loci and annotation preparation → <strong>STOP 2: candidate models, alignment and tree</strong> → STOP 3: synteny / genomic context → STOP 4: final dossier.</p><p>Annotation or Tiberius models are prepared before selection when needed. Choose several models per locus. Rebuild the alignment and tree together, then inspect synteny. Initial checkboxes are comparison suggestions, never approved models. Existing manual decisions remain separate.</p><p>Protein and CDS views remain distinct. Missing query or annotation CDS is not invented. Similar isoforms from one locus are alternative models, not independent species observations.</p></main></body></html>'
    html = attach(html, pool)
    html = html.replace('</body>', '<script>'+Path(__file__).with_name('alignment_review.txt').read_text()+'</script></body>')
    if legacy and provenance.get('checkpoint'):
        html = prepared_page(Path(provenance['checkpoint']), pool, title)
    (root/'report.html').write_text(html)
    dump(root/'manifest.json', dict(schema=SCHEMA, status='awaiting_model_review', created=now(),
        review_id=pool['review_id'], source=provenance, models=len([s for s in pool['sequences'] if not s['reference']]),
        references=len([s for s in pool['sequences'] if s['reference']]),
        pool_sha256=digest([[s['id'], s['protein'], s['cds'], s['group']] for s in pool['sequences']]),
        biological_approval_asserted=False, tiberius_run=False))
    files = ['review_data.json','report.html','manifest.json']
    if (root/'review_app.zip').is_file():
        files.append('review_app.zip')
    dump(root/'seal.json', {name:sha256(root/name) for name in files})
    return root, pool


def validate_selection(pool, selection):
    if selection.get('schema') != 'alignment_model_selection_v1' or selection.get('review_id') != pool['review_id']:
        raise ValueError('Wrong alignment-model review')
    if selection.get('source_review_id') != pool['source_review_id']:
        raise ValueError('Wrong source checkpoint')
    ids = selection.get('ids')
    by_id = {s['id']:s for s in pool['sequences']}
    if not isinstance(ids, list) or not ids or len(set(ids)) != len(ids) or set(ids)-set(by_id):
        raise ValueError('Unknown, duplicate or empty model set')
    if not any(not by_id[sid]['reference'] for sid in ids):
        raise ValueError('No candidate model selected')
    if selection.get('biological_approval_asserted') is not False:
        raise ValueError('Alignment membership cannot assert biological approval')
    return [by_id[sid] for sid in ids]
