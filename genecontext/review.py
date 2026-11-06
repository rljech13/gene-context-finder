"""Explicit, file-backed review gates. Blank decisions never approve a candidate."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from html import escape
import json
from pathlib import Path
import shutil
import sys
from uuid import uuid4

from .comparison import (make_groups, orient_interval, same_locus, select_candidates, table)
from .sequence import Model, coding_sequence, read_fasta, read_gff, reverse_complement, sha256, translate, union_length, write_fasta, write_tsv
from .protein_discovery import SCORE_FIELDS
from .tools import Runner, tool_versions


def optional_fasta(path):
    return read_fasta(path) if path.exists() and path.stat().st_size else {}


def now():
    return datetime.now(timezone.utc).isoformat()


def new_root(output, source):
    root, source = Path(output).resolve(), Path(source).resolve()
    if root == source or source in root.parents or root in source.parents:
        raise ValueError('Review output must be separate from its source')
    if root.exists() and any(root.iterdir()):
        raise ValueError('Review output must be new or empty')
    root.mkdir(parents=True, exist_ok=True)
    return root


def seal(root):
    # Large scratch databases/indexes are not review artifacts and are not copied.
    directories = {'windows', 'context_prediction', 'context_predictions', 'evidence', 'figures'}
    files = [f for f in root.iterdir() if f.is_file() and f.name != 'seal.json'
             and not f.name.endswith('.decisions.tsv')]
    files += [f for name in directories for f in (root/name).rglob('*')
              if f.is_file() and 'cache' not in f.relative_to(root/name).parts]
    (root/'seal.json').write_text(json.dumps({str(f.relative_to(root)): sha256(f) for f in sorted(files)}, indent=2)+'\n')


def verify(root, stage, source_override=None):
    root = Path(root).resolve()
    hashes = json.loads((root/'seal.json').read_text())
    for name, digest in hashes.items():
        path = (root/name).resolve()
        if root not in path.parents or not path.is_file() or sha256(path) != digest:
            raise ValueError('Review evidence changed or is missing: '+name)
    manifest = json.loads((root/'manifest.json').read_text())
    if manifest.get('stage') != stage or manifest.get('status') != 'awaiting_review':
        raise ValueError('Expected an awaiting-review '+stage+' checkpoint')
    if not manifest.get('review_id'):raise ValueError('Checkpoint lacks a review ID; regenerate it from saved results')
    source = Path(source_override).resolve() if source_override else Path(manifest['source_results'])
    if source_override:manifest.update(recorded_source_results=manifest['source_results'],source_results=str(source))
    for name, digest in manifest['source_checksums'].items():
        if not (source/name).is_file() or sha256(source/name) != digest:
            raise ValueError('Source search evidence changed: '+name)
    return root, manifest


def source_checksums(source):
    files = [source/n for n in ['manifest.json', 'models.tsv', 'dna_regions.tsv']]
    files += [p for p in (source/'sequences').rglob('*') if p.name in {'cds.fna', 'protein.faa', 'dna_regions.fna'}]
    return {str(p.relative_to(source)): sha256(p) for p in files if p.is_file()}


def decisions(path, keys, expected, actions, review_id=None):
    rows = table(path)
    result = {}
    for row in rows:
        if review_id is not None and row.get('review_id','').strip()!=review_id:
            raise ValueError('Decision file belongs to another review checkpoint')
        if not all(k in row for k in [*keys, 'decision']):
            raise ValueError('Decision TSV requires '+', '.join([*keys, 'decision']))
        key = tuple(row[k].strip() for k in keys)
        if key not in expected or key in result:
            raise ValueError('Unknown or duplicate decision key: '+str(key))
        row = {k: (v or '').strip() for k, v in row.items()}
        if row['decision'] not in actions:
            raise ValueError('Pending/invalid decision for '+str(key)+'; choose '+', '.join(sorted(actions)))
        if row['decision'] in {'exclude', 'unresolved'} and not row.get('notes'):
            raise ValueError('Exclusions and unresolved cases require notes: '+str(key))
        result[key] = row
    if set(result) != set(expected):
        raise ValueError('Missing explicit decisions for '+str(sorted(set(expected)-set(result))))
    return result


def page(path, title, text, rows, columns):
    cells = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>'+escape(title)+'</title>',
             '<style>body{font:15px system-ui;margin:25px}table{border-collapse:collapse}th,td{border:1px solid #ccd;padding:7px;text-align:left;vertical-align:top}code{background:#eef2f6}aside{padding:15px;background:#fff7ed}td{max-width:600px;overflow-wrap:anywhere}</style>',
             '<h1>'+escape(title)+'</h1>', text, '<table><tr>'+''.join('<th>'+escape(c)+'</th>' for c in columns)+'</tr>']
    for row in rows:
        cells.append('<tr>'+''.join('<td>'+escape(str(row.get(c, '')))+'</td>' for c in columns)+'</tr>')
    cells.append('</table></html>');path.write_text('\n'.join(cells))


def review_search(args):
    from . import __version__
    source = Path(args.results).resolve()
    original = json.loads((source/'manifest.json').read_text())
    groups = make_groups(original['queries'], args.query_groups)
    _, candidates, observations = select_candidates(groups, original['genomes'], table(source/'models.tsv'), table(source/'dna_regions.tsv'))
    root = new_root(args.out, source)
    for query in original['queries']:
        for name in ['protein.faa', 'cds.fna', 'dna_regions.fna']:
            src = source/'sequences'/query['id']/name
            if src.exists():
                dst = root/'evidence'/query['id']/name;dst.parent.mkdir(parents=True, exist_ok=True);shutil.copyfile(src, dst)
    source_rows={r['model_id']:r for r in table(source/'models.tsv')}
    group_names={g['id']:g['name'] for g in groups}
    for c in candidates:
        c['suggested'] = c.pop('selected')
        c['group_name']=group_names[c['group']]
        c.update({k:source_rows.get(c['model_id'],{}).get(k,'') for k in ['protein_aa','cds_bp','issues','alignment_query_coverage']})
    columns = ['group', 'group_name', 'genome_id', 'model_id', 'locus_rank', 'representative', 'suggested', 'query', 'contig', 'start', 'end', 'strand', 'status', 'protein_aa', 'cds_bp', 'issues', 'source', 'source_gene', 'source_model', 'coverage', 'identity', 'alignment_query_coverage', 'ranking_score']
    write_tsv(root/'candidates.tsv', candidates, columns)
    review_id=str(uuid4())
    template = [dict(review_id=review_id,group=o['group'], genome_id=o['genome_id'], decision='', model_ids='', notes='',
                     suggested_model_id=o['model_id'], search_status=o['status'], candidate_loci=o['candidate_loci']) for o in observations]
    write_tsv(root/'loci.decisions.tsv', template, ['review_id', 'group', 'genome_id', 'decision', 'model_ids', 'notes', 'suggested_model_id', 'search_status', 'candidate_loci'])
    manifest = dict(version=__version__, review_id=review_id, stage='locus_review', status='awaiting_review', created=now(),
                    source_results=str(source), source_checksums=source_checksums(source), groups=groups,
                    pairs=[dict(group=o['group'], genome_id=o['genome_id']) for o in observations],
                    query_groups_sha256=sha256(args.query_groups) if args.query_groups else None)
    (root/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    text = '<aside><strong>STOP 1: no locus is approved.</strong> Rankings are suggestions, not gene-family or copy assignments.</aside><p>Edit <a href="loci.decisions.tsv">loci.decisions.tsv</a>: <code>inspect</code> with one or more exact model IDs separated by semicolons; <code>exclude</code> or <code>unresolved</code> with notes. Every gene/genome pair needs an explicit decision, including pairs without candidates. Several loci can proceed to context inspection. Blank rows block the next stage.</p><p><a href="candidates.tsv">All candidates and isoforms</a> · <a href="manifest.json">Source fingerprints and groups</a>. FASTA evidence is in evidence/&lt;query_id&gt;/.</p>'
    if original.get('errors'):text+='<h2>Search errors</h2><pre>'+escape(json.dumps(original['errors'],indent=2))+'</pre>'
    page(root/'report.html', 'Stage 1 — review candidate loci', text, candidates, columns)
    seal(root);print('STOP 1 — review '+str(root/'report.html'), file=sys.stderr)
    return 0


def search(args):
    from .pipeline import run
    if args.word_size < 4 or args.protein_evalue <= 0 or args.dna_evalue <= 0:
        raise ValueError('Word size must be >=4 and E-values must be positive')
    output = Path(args.out).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Search output must be new or empty')
    output.mkdir(parents=True, exist_ok=True)
    searching = argparse.Namespace(**vars(args));searching.out = str(output/'search')
    code = run(searching)
    reviewing = argparse.Namespace(results=searching.out, out=str(output/'review'), query_groups=args.query_groups)
    review_search(reviewing)
    (output/'report.html').write_text('<!doctype html><meta charset="utf-8"><h1>Search completed — manual review required</h1><p><a href="review/report.html">STOP 1: candidate review and decisions</a> · <a href="search/report.html">Search evidence</a></p>')
    return code


def contexts(args):
    from .comparison import compare
    root, checkpoint = verify(args.review, 'locus_review',getattr(args,'search_results',None))
    source = Path(checkpoint['source_results']);original = json.loads((source/'manifest.json').read_text())
    groups = checkpoint['groups'];expected = {(r['group'], r['genome_id']) for r in checkpoint['pairs']}
    choices = decisions(args.decisions, ['group', 'genome_id'], expected, {'inspect', 'exclude', 'unresolved'}, checkpoint['review_id'])
    selected, candidates, observations = select_candidates(groups, original['genomes'], table(source/'models.tsv'), table(source/'dna_regions.tsv'))
    # Rebuild the selection from explicit IDs. The ranked default never advances.
    all_rows = {r['model_id']: r for r in table(source/'models.tsv')}
    selected = [];allowed = {(r['group'], r['genome_id'], r['model_id']) for r in candidates}
    for o in observations:
        key=(o['group'], o['genome_id']);choice=choices[key];ids=[x.strip() for x in choice.get('model_ids', '').split(';') if x.strip()]
        if choice['decision'] != 'inspect':
            if ids:raise ValueError('Only inspect decisions may contain model_ids: '+str(key))
            o.update(status='review_'+choice['decision'], model_id='', error=choice['notes']);continue
        if not ids or len(ids)!=len(set(ids)):raise ValueError('inspect requires distinct model_ids: '+str(key))
        for mid in ids:
            if (*key, mid) not in allowed:raise ValueError('Model does not belong to decision pair: '+mid)
            row=all_rows.get(mid)
            if row is None:
                reg=next(r for r in table(source/'dna_regions.tsv') if r['region_id']==mid)
                row=dict(reg, model_id=mid, source_model=mid, source_gene='', status='dna_region', cds_exons='[]', protein_aa=0, source=reg.get('source') or 'blastn')
            if any(r['group']==key[0] and r['genome_id']==key[1] and same_locus(row,r) for r in selected):
                raise ValueError('Inspect one seed per locus; isoforms are available at STOP 2: '+mid)
            selected.append(dict(row, group=key[0], locus_id=f'L{len(selected)+1:06d}', selection_status='manual_inspection', candidate_loci=o['candidate_loci']))
        o.update(status='manual_inspection', model_id=';'.join(ids))
    if not selected:return terminal_review(args.out,root,'locus_review',choices,args.decisions)
    selected_ids={r['model_id'] for r in selected}
    for c in candidates:c['selected']=c['model_id'] in selected_ids
    if args.query_groups or args.select_loci:raise ValueError('Grouping and locus selection are fixed by STOP 1 decisions')
    # Preserve the exact reviewed grouping without a temporary mutable mapping file.
    args.results=str(source);args._review_groups=groups
    args._locus_review=dict(review_id=checkpoint['review_id'],path=str(root), decisions=str(Path(args.decisions).resolve()), decisions_sha256=sha256(args.decisions), choices=list(choices.values()))
    return compare(args, selection=(selected, candidates, observations), context_only=True)


def terminal_review(output, source, stage, choices, decision_file):
    from . import __version__
    root=new_root(output,source)
    shutil.copyfile(decision_file,root/'decisions.tsv')
    rows=list(choices.values())
    manifest=dict(version=__version__,stage=stage,status='completed_no_comparable_targets',finished=now(),review_source=str(source),decisions_sha256=sha256(decision_file),choices=rows,commands=[])
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    page(root/'report.html','Review complete — no targets retained','<p>Every pair was explicitly excluded or left unresolved. No search, prediction or alignment was run. These decisions do not prove biological gene absence. <a href="decisions.tsv">Recorded decisions</a></p>',rows,['group','genome_id','decision','notes'])
    return 0


def _copy_evidence(source, root):
    for name in ['windows', 'context_prediction', 'context_predictions']:
        if (source/name).exists():
            shutil.copytree(source/name, root/name, ignore=shutil.ignore_patterns('cache'))
    for original,saved in [('target_review.tsv','prior_target_identity_review.tsv'),('identity_audit.tsv','prior_protein_identity_audit.tsv')]:
        prior=source/'validation'/original
        if prior.exists():shutil.copyfile(prior,root/saved)
        elif (source/saved).exists():shutil.copyfile(source/saved,root/saved)
    for name in ['provided_neighbors.tsv', 'context_prediction.json', 'context_reference.json', 'context_predictions.tsv', 'context_prediction_cds.fna', 'context_prediction_proteins.faa']:
        if (source/name).exists():shutil.copyfile(source/name, root/name)


def _prediction_features(root, selected, context, proteins, metadata):
    """Keep all predicted models reviewable, including the old target-overlap models."""
    path=root/'context_predictions.tsv'
    if not path.exists():return context
    rows={r['locus_id']:r for r in selected};aa=optional_fasta(root/'context_prediction_proteins.faa')
    result=[]
    for r in table(path):
        row=rows[r['locus_id']];pid=r['protein_id'];start,end=int(r['start']),int(r['end'])
        usable=r['coding_status']!='disrupted' and len(aa.get(pid,''))>=30 and '*' not in aa.get(pid,'')
        if usable:
            proteins[pid]=aa[pid];metadata.setdefault(pid,dict(genome_id=row['genome_id'],gene_id=pid,name=r.get('name') or 'uncharacterized',evidence='ab_initio_prediction',coding_status=r['coding_status']))
        a,b=orient_interval(row,start,end)
        result.append(dict(locus_id=row['locus_id'],group=row['group'],genome_id=row['genome_id'],gene_id=pid,protein_id=pid if usable else '',name=r.get('name') or 'uncharacterized',contig=row['contig'],start=start,end=end,strand=r['strand'],relative_start=a,relative_end=b,relative_strand='+' if r['strand']==row['strand'] else '-',evidence='ab_initio_prediction',clipped=r['coding_status']=='partial',coding_status=r['coding_status'],identification_status=r.get('identification_status') or 'not_identified',reference_protein=r.get('reference_protein',''),reference_identity=r.get('reference_identity',''),reference_coverage=r.get('reference_coverage','')))
    return result


def _options(source, root, original, selected):
    source_models=table(source/'models.tsv');source_seqs={}
    for q in original['queries']:
        for kind,file in [('cds','cds.fna'),('protein','protein.faa')]:
            path=source/'sequences'/q['id']/file
            source_seqs[q['id'],kind]=optional_fasta(path)
    preds={}
    for path in (root/'context_predictions').glob('*.gff3'):
        for model in read_gff(path)[0]:preds[model.id.removesuffix('.t1')]=model
    raw=table(root/'context_predictions.tsv') if (root/'context_predictions.tsv').exists() else []
    pred_cds=read_fasta(root/'context_prediction_cds.fna') if raw else {};pred_aa=optional_fasta(root/'context_prediction_proteins.faa') if raw else {}
    options=[];sequences={}
    for locus in selected:
        if locus.get('comparison_error') or 'window_start' not in locus:continue
        lid=locus['locus_id']
        for row in source_models:
            if row['genome_id']!=locus['genome_id'] or row['query']!=locus['query'] or not same_locus(row,locus):continue
            parts=json.loads(row['cds_exons']);start,end=int(row['start']),int(row['end'])
            if start<locus['window_start'] or end>locus['window_end']:continue
            oid=lid+'.source.'+row['model_id'];cds=source_seqs[row['query'],'cds'][row['model_id']];aa=source_seqs[row['query'],'protein'][row['model_id']]
            item=dict(option_id=oid,locus_id=lid,group=locus['group'],genome_id=locus['genome_id'],origin='source',model_id=row['model_id'],query=row['query'],start=start,end=end,strand=row['strand'],status=row['status'],issues=row.get('issues',''),protein_aa=len(aa),native_cds=parts,eligible_target=row['status']!='disrupted',name=locus['target_name'])
            item.update(source_model=row.get('source_model',''),source_gene=row.get('source_gene',''),scores={k:row.get(k,'') for k in ['coverage','identity',*SCORE_FIELDS]})
            ordered=sorted(parts,reverse=row['strand']=='-');extra=len(cds)-sum(b-a for a,b,*_ in parts)+(ordered[0][2] or 0)
            terminal=min(a for a,b,*_ in parts) if row['strand']=='-' else max(b for a,b,*_ in parts)
            item['native_stops']=[(terminal-extra,terminal) if row['strand']=='-' else (terminal,terminal+extra)] if 0<extra<=3 else []
            options.append(item);sequences[oid]=dict(cds=cds,protein=aa)
        for rec in raw:
            if rec['locus_id']!=lid:continue
            pid=rec['protein_id'];model=preds[pid];oid=lid+'.prediction.'+pid
            item=dict(option_id=oid,locus_id=lid,group=locus['group'],genome_id=locus['genome_id'],origin='prediction',model_id=pid,query=locus['query'],start=int(rec['start']),end=int(rec['end']),strand=rec['strand'],status=rec['coding_status'],issues=rec['issues'],protein_aa=len(pred_aa[pid]),native_cds=model.cds,native_stops=model.stops,eligible_target=rec['strand']==locus['strand'] and rec['coding_status']!='disrupted',name=rec.get('name') or 'uncharacterized')
            options.append(item);sequences[oid]=dict(cds=pred_cds[pid],protein=pred_aa[pid])
        options.append(dict(option_id=lid+'.dna_only',locus_id=lid,group=locus['group'],genome_id=locus['genome_id'],origin='dna_only',model_id=locus['model_id'],query=locus['query'],start=locus['target_start'],end=locus['target_end'],strand=locus['strand'],status='dna_region',issues='',protein_aa=0,native_cds=[],eligible_target=True,name=locus['target_name']))
    return options,sequences


def _exon_figure(path, locus, options):
    width=1400;scale=1200/locus['window_bp'];parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="HEIGHT" viewBox="0 0 1400 HEIGHT"><rect width="100%" height="100%" fill="white"/>', '<style>text{font:12px sans-serif;fill:#172033}</style>',f'<text x="20" y="24">{escape(locus["genome_id"])} | {escape(locus["locus_id"])} | {escape(locus["contig"])}</text>',f'<text x="20" y="46">Native window {locus["window_start"]}–{locus["window_end"]} (0-based, half-open). Rows oriented with discovery target. Boxes = CDS.</text>'];y=80
    for op in options:
        if op['origin']=='dna_only':continue
        color='#2563eb' if op['origin']=='source' else '#059669'
        intervals=[orient_interval(locus,a,b) for a,b,*_ in op['native_cds']]
        if intervals:
            a=min(a for a,b in intervals);b=max(b for a,b in intervals)
            parts.append(f'<path d="M{100+a*scale:.2f} {y}H{100+b*scale:.2f}" stroke="{color}"/>')
            for a,b in intervals:parts.append(f'<rect x="{100+a*scale:.2f}" y="{y-5}" width="{max(1,(b-a)*scale):.2f}" height="10" fill="{color}"/>')
        parts.append(f'<text x="20" y="{y+20}">{escape(op["option_id"])} | {op["protein_aa"]} aa | {escape(op["status"])} | strand {op["strand"]}</text>');y+=54
    parts.append('</svg>');path.write_text('\n'.join(parts).replace('HEIGHT',str(y+20)))


def save_context_review(args, source, root, original, groups, selected, candidates, observations, context, proteins, metadata, windows, runner, errors, manifest):
    # No DNA/target-protein alignments are run before this second gate.
    options,sequences=_options(source,root,original,selected)
    by_lid={r['locus_id']:r for r in selected}
    prior=table(root/'prior_protein_identity_audit.tsv') if (root/'prior_protein_identity_audit.tsv').exists() else []
    prior={r['model_id']:r for r in prior if r.get('rank')=='1'}
    group_names={g['id']:g['name'] for g in groups}
    for op in options:
        hit=prior.get(op['model_id'],{})
        op.update(prior_reference_gene=hit.get('gene',''),prior_candidate_coverage=hit.get('query_coverage',''),prior_reference_coverage=hit.get('reference_coverage',''))
        op['group_name']=group_names[op['group']]
        loc=by_lid[op['locus_id']]
        a,b=(loc['target_start'],loc['target_end']) if op['origin']=='source' else (op['start'],op['end'])
        op.update(left_flank_bp=a-loc['window_start'],right_flank_bp=loc['window_end']-b)
        op['full_requested_flanks']=min(op['left_flank_bp'],op['right_flank_bp'])>=args.window
    all_context=_prediction_features(root,selected,context,proteins,metadata)
    from .context_compare import homology_groups
    from .comparison_report import synteny_svg
    if getattr(args,'_imported_contexts',False):
        families={f['protein_id']:f.get('homology_group','') for f in context}
        for f in all_context:f['homology_group']=families.get(f['protein_id'],'')
    else:
        try:
            families,_=homology_groups(args,root,all_context,proteins,metadata,runner)
            for f in all_context:f['homology_group']=families.get(f['protein_id'],'')
        except (OSError,ValueError,RuntimeError) as error:errors.append(dict(stage='provisional_neighbor_homology',error=str(error)))
    (root/'figures').mkdir()
    for row in selected:
        if not row.get('comparison_error'):_exon_figure(root/'figures'/(row['locus_id']+'.svg'),row,[o for o in options if o['locus_id']==row['locus_id']])
    for group in groups:synteny_svg(root/'figures'/(group['id']+'.svg'),group['name']+' — provisional, unreviewed contexts',[r for r in selected if r['group']==group['id'] and not r.get('comparison_error')],all_context)
    columns=['option_id','locus_id','group','group_name','genome_id','origin','model_id','status','issues','protein_aa','strand','start','end','eligible_target','left_flank_bp','right_flank_bp','full_requested_flanks','name','prior_reference_gene','prior_candidate_coverage','prior_reference_coverage']
    write_tsv(root/'model_options.tsv',options,columns)
    write_fasta(root/'model_options.cds.fna',[(k,v['cds']) for k,v in sequences.items()]);write_fasta(root/'model_options.proteins.faa',[(k,v['protein']) for k,v in sequences.items()])
    review_id=str(uuid4())
    rows=[dict(review_id=review_id,group=o['group'],genome_id=o['genome_id'],decision='',option_id='',target_name='',allow_partial='',allow_short_flanks='',notes='') for o in observations]
    write_tsv(root/'models.decisions.tsv',rows,['review_id','group','genome_id','decision','option_id','target_name','allow_partial','allow_short_flanks','notes'])
    nr=[dict(review_id=review_id,locus_id=f['locus_id'],gene_id=f['gene_id'],decision='',name='',notes='',suggested_name=f['name']) for f in all_context]
    write_tsv(root/'neighbors.decisions.tsv',nr,['review_id','locus_id','gene_id','decision','name','notes','suggested_name'])
    write_tsv(root/'neighbor_options.tsv',all_context,['locus_id','gene_id','name','start','end','strand','coding_status','identification_status','reference_identity','reference_coverage'])
    state=dict(original=original,groups=groups,selected=selected,candidates=candidates,observations=observations,context=all_context,proteins=proteins,metadata=metadata,windows=windows,options=options,sequences=sequences,parameters=vars(args),errors=errors)
    (root/'state.json').write_text(json.dumps(state,indent=2,default=str)+'\n')
    from . import __version__
    manifest.update(version=__version__,review_id=review_id,review_source_code_sha256={p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')},stage='context_review',status='awaiting_review',finished=now(),commands=runner.commands,errors=errors,
                    source_checksums=source_checksums(source),locus_review=getattr(args,'_locus_review',None))
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2,default=str)+'\n')
    links=' · '.join('<a href="figures/'+escape(g['id'])+'.svg">'+escape(g['name'])+' contexts</a>' for g in groups)
    text='<aside><strong>STOP 2: targets and neighbours are unapproved.</strong> No target protein/CDS or genomic DNA alignments have been run by this stage.</aside><p>Edit <a href="models.decisions.tsv">models.decisions.tsv</a>: <code>keep</code> with exactly one option ID per gene/genome pair; or <code>exclude</code>/<code>unresolved</code> with notes. A partial CDS requires <code>allow_partial=yes</code> and notes. Short flanks require <code>allow_short_flanks=yes</code> and notes. A <code>dna_only</code> option permits DNA/context comparison without a fabricated protein.</p><p>Edit <a href="neighbors.decisions.tsv">neighbors.decisions.tsv</a>: <code>keep</code>, <code>unnamed</code>, or <code>exclude</code>. Only neighbours of retained targets require decisions; models overlapping the chosen target CDS are automatically treated as target models. Renaming or exclusion requires notes.</p><p><a href="model_options.tsv">Target model options</a> · <a href="model_options.cds.fna">CDS</a> · <a href="model_options.proteins.faa">Proteins</a> · <a href="neighbor_options.tsv">Neighbour evidence</a> · '+links+'</p><p>Per-locus exon tracks: '+ ' · '.join('<a href="figures/'+r['locus_id']+'.svg">'+r['locus_id']+'</a>' for r in selected if not r.get('comparison_error'))+'</p><p>Choosing a new model preserves the prepared DNA window. Actual flank sizes and any shortage are recorded in the final report. Predictions are hypotheses; selecting a row records a human decision, not automatic proof of orthology.</p>'
    if getattr(args,'_imported_contexts',False):text+='<aside>Imported automatic selections remain unapproved. To inspect alternative loci absent from these saved windows, complete STOP 1 and run contexts for those choices.</aside>'
    if (root/'prior_target_identity_review.tsv').exists():text+='<p><a href="prior_target_identity_review.tsv">Previous reference-identity flags for the original automatic targets</a>. These flags do not classify newly chosen models.</p>'
    if (root/'prior_protein_identity_audit.tsv').exists():text+='<p><a href="prior_protein_identity_audit.tsv">Previous independent BLASTP evidence for source and predicted proteins</a>. Reference-name/coverage columns are supporting evidence, not orthology assignments.</p>'
    if errors:text+='<h2>Processing errors</h2><pre>'+escape(json.dumps(errors,indent=2))+'</pre>'
    page(root/'report.html','Stage 2 — review target CDS and neighbours',text,options,columns)
    seal(root);print('STOP 2 — review '+str(root/'report.html'),file=sys.stderr)
    return 1 if errors or original.get('errors') else 0


def _target_neighbor(feature, option, row):
    if feature['strand']!=option['strand']:return False
    if not option['native_cds']:
        return feature['gene_id']==row.get('source_gene')
    from .context_compare import _overlap
    # Raw predictions have exact CDS; provided reference features use their envelope.
    parts=feature.get('_native_cds') or [(feature['start'],feature['end'])]
    other=[(a,b) for a,b,*_ in option['native_cds']]
    return _overlap(parts,other)>=.8*min(union_length(parts),union_length(other))


def validate_model_sequence(option, row, window, expected):
    class WindowGenome:
        def fetch(self, contig, start, end, strand='+'):
            if contig!=row['contig'] or not row['window_start']<=start<=end<=row['window_end']:
                raise ValueError('Chosen model extends beyond the prepared window')
            a,b=orient_interval(row,start,end);piece=window[a:b].upper()
            return piece if strand==row['strand'] else reverse_complement(piece)
    model=Model(option['model_id'],row['contig'],option['strand'],cds=[tuple(p) for p in option['native_cds']],stops=[tuple(p) for p in option.get('native_stops',[])])
    cds,protein,_,_=coding_sequence(model,WindowGenome())
    if cds!=expected['cds'] or protein!=expected['protein']:
        raise ValueError('Chosen model does not round-trip to its window DNA: '+option['option_id'])


def write_approved_gff(root, selected, options):
    from urllib.parse import quote
    with open(root/'reviewed_targets.gff3','w') as handle:
        handle.write('##gff-version 3\n')
        for row in selected:
            op=options[row['locus_id']];pid=row['locus_id'];common=row['contig']+'\tmanual_review\t'
            kind='gene' if op['native_cds'] else 'region'
            handle.write(common+f'{kind}\t{op["start"]+1}\t{op["end"]}\t.\t{op["strand"]}\t.\tID={pid};Name={quote(row["target_name"],safe="")}\n')
            if not op['native_cds']:continue
            handle.write(common+f'mRNA\t{op["start"]+1}\t{op["end"]}\t.\t{op["strand"]}\t.\tID={pid}.t1;Parent={pid}\n')
            for a,b,phase in op['native_cds']:
                handle.write(common+f'CDS\t{a+1}\t{b}\t.\t{op["strand"]}\t{phase if phase is not None else "."}\tParent={pid}.t1\n')
            for a,b in op.get('native_stops',[]):
                handle.write(common+f'stop_codon\t{a+1}\t{b}\t.\t{op["strand"]}\t0\tParent={pid}.t1\n')


def compare_reviewed(args):
    from . import __version__
    from .comparison import finish_comparison
    checkpoint,manifest=verify(args.review,'context_review',getattr(args,'search_results',None))
    state=json.loads((checkpoint/'state.json').read_text());source=Path(manifest['source_results'])
    expected={(r['group'],r['genome_id']) for r in state['observations']}
    choices=decisions(args.models,['group','genome_id'],expected,{'keep','exclude','unresolved'},manifest['review_id'])
    options={o['option_id']:o for o in state['options']};loci={r['locus_id']:r for r in state['selected']}
    selected=[];sequence_overrides={};selected_options={}
    for key,choice in choices.items():
        oid=choice.get('option_id','')
        if choice['decision']!='keep':
            if oid:raise ValueError('Only keep decisions may contain option_id: '+str(key))
            continue
        if oid not in options:raise ValueError('Unknown option_id: '+oid)
        op=options[oid]
        if (op['group'],op['genome_id'])!=key or not op['eligible_target']:raise ValueError('Ineligible/wrong-pair target model: '+oid)
        row=dict(loci[op['locus_id']])
        if not op['full_requested_flanks'] and (choice.get('allow_short_flanks')!='yes' or not choice.get('notes')):raise ValueError('Short flanks require allow_short_flanks=yes and notes: '+oid)
        if row.get('comparison_error'):raise ValueError('Cannot keep a failed window: '+oid)
        if op['status']=='partial' and (choice.get('allow_partial')!='yes' or not choice.get('notes')):raise ValueError('Partial CDS requires allow_partial=yes and notes: '+oid)
        row.update(discovery_model_id=row['model_id'],model_id=op['model_id'],protein_aa=op['protein_aa'],status=op['status'],target_name=choice.get('target_name') or op['name'],selection_status='manually_reviewed',review_option_id=oid,review_origin=op['origin'])
        row.update(op.get('scores') or {k:'' for k in ['coverage','identity',*SCORE_FIELDS]})
        row['source_model']=op.get('source_model') or op['model_id']
        if op['origin']!='source':row.update(target_start=op['start'],target_end=op['end'],anchor_kind='reviewed_'+op['origin'],source='reviewed_'+op['origin'],source_gene=op['model_id'])
        row['relative_cds']=[orient_interval(row,a,b) for a,b,*_ in op['native_cds']]
        row['target_relative_start'],row['target_relative_end']=orient_interval(row,row['target_start'],row['target_end'])
        radius=state['parameters']['window'];row.update(left_clipped=row['window_start']>row['target_start']-radius,right_clipped=row['window_end']<row['target_end']+radius,
            actual_left_flank_bp=row['target_start']-row['window_start'],actual_right_flank_bp=row['window_end']-row['target_end'])
        if op['origin']!='dna_only':
            seq=state['sequences'][oid]
            validate_model_sequence(op,row,state['windows'][row['locus_id']],seq)
            if len(seq['cds'])%3 or translate(seq['cds']).removesuffix('*')!=seq['protein']:raise ValueError('Chosen CDS/protein inconsistency: '+oid)
            sequence_overrides[row['locus_id']]=seq
        selected.append(row);selected_options[row['locus_id']]=op
    if not selected:return terminal_review(args.out,checkpoint,'context_review',choices,args.models)
    # Recompute target membership from exact predicted CDS, not the old target flag.
    predicted={}
    for path in (checkpoint/'context_predictions').glob('*.gff3'):
        for model in read_gff(path)[0]:predicted[model.id.removesuffix('.t1')]=[(a,b) for a,b,*_ in model.cds]
    needed=[]
    for f in state['context']:
        lid=f['locus_id']
        if lid not in selected_options:continue
        f['_native_cds']=predicted.get(f['gene_id'])
        if not _target_neighbor(f,selected_options[lid],loci[lid]):needed.append(f)
    all_neighbor_keys={(f['locus_id'],f['gene_id']) for f in state['context']}
    # Unchosen windows and target-overlap rows may stay blank, but unknown keys never pass.
    raw=table(args.neighbors);nrows={}
    for r in raw:
        if r.get('review_id','').strip()!=manifest['review_id']:raise ValueError('Neighbour decisions belong to another review checkpoint')
        key=(r.get('locus_id','').strip(),r.get('gene_id','').strip())
        if key not in all_neighbor_keys or key in nrows:raise ValueError('Unknown/duplicate neighbour decision: '+str(key))
        nrows[key]={k:(v or '').strip() for k,v in r.items()}
    context=[]
    for f in needed:
        key=(f['locus_id'],f['gene_id']);choice=nrows.get(key,{})
        if choice.get('decision') not in {'keep','unnamed','exclude'}:raise ValueError('Pending neighbour decision: '+str(key))
        if (choice['decision']=='exclude' or choice.get('name')) and not choice.get('notes'):raise ValueError('Neighbour exclusion/rename requires notes: '+str(key))
        if choice['decision']=='exclude':continue
        f['name']='uncharacterized' if choice['decision']=='unnamed' else choice.get('name') or f['name']
        f['identification_status']='manually_reviewed_unnamed' if choice['decision']=='unnamed' else 'manually_reviewed'
        f.pop('_native_cds',None);context.append(f)
    params=argparse.Namespace(**state['parameters']);params.out=args.out;params.threads=args.threads
    for name in ['dna_engine','dna_min_block','dna_min_identity','reference_genome','max_alignment_bases']:
        if getattr(args,name,None) is not None:setattr(params,name,getattr(args,name))
    params.include_partial=any(r['status']=='partial' for r in selected)
    versions=tool_versions(['blastp','makeblastdb','mafft',params.dna_engine])
    if any(not v['path'] for v in versions.values()):raise ValueError('Missing reviewed-comparison tools')
    if params.reference_genome and params.reference_genome not in {g['genome_id'] for g in state['original']['genomes']}:raise ValueError('Unknown --reference-genome')
    root=new_root(args.out,checkpoint);_copy_evidence(checkpoint,root)
    proteins=dict(state['proteins']);metadata=dict(state['metadata'])
    live={f['protein_id'] for f in context if f['protein_id']}
    supporting={metadata[p].get('reference_protein') for p in live if p in metadata}
    proteins={p:s for p,s in proteins.items() if p in live or p in supporting};metadata={p:m for p,m in metadata.items() if p in proteins}
    for f in context:
        if f['protein_id'] in metadata:metadata[f['protein_id']]['name']=f['name']
    observations=state['observations']
    for o in observations:
        choice=choices[o['group'],o['genome_id']];o.update(status='review_'+choice['decision'],model_id='',error=choice.get('notes',''))
        row=next((r for r in selected if (r['group'],r['genome_id'])==(o['group'],o['genome_id'])),None)
        if row:o.update(model_id=row['model_id'])
    chosen_ids={r['model_id'] for r in selected}
    for c in state['candidates']:c['selected']=c['model_id'] in chosen_ids
    runner=Runner(root);errors=list(state['errors'])
    output_manifest=dict(version=__version__,stage='reviewed_comparison',started=now(),status='running',source_results=str(source),groups=state['groups'],parameters=vars(params),tools=versions,
        source_code_sha256={p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')},
        review=dict(context_checkpoint=str(checkpoint),checkpoint_seal_sha256=sha256(checkpoint/'seal.json'),locus_review=manifest.get('locus_review'),model_decisions_sha256=sha256(args.models),neighbor_decisions_sha256=sha256(args.neighbors),prediction_reused=True,model_decisions=list(choices.values())),context_prediction=manifest.get('context_prediction',{}))
    for src,name in [(args.models,'approved_models.tsv'),(args.neighbors,'approved_neighbors.tsv')]:shutil.copyfile(src,root/name)
    prediction_table=root/'context_predictions.tsv'
    if prediction_table.exists():
        shutil.copyfile(prediction_table,root/'original_context_predictions.tsv')
        raw_predictions=table(prediction_table);cols=[*raw_predictions[0], 'original_target_overlap'] if raw_predictions else []
        for rec in raw_predictions:
            rec['original_target_overlap']=rec['target_overlap']
            lid=rec['locus_id']
            if lid in selected_options:
                f=dict(gene_id=rec['protein_id'],strand=rec['strand'],start=int(rec['start']),end=int(rec['end']),_native_cds=predicted.get(rec['protein_id']))
                rec['target_overlap']=_target_neighbor(f,selected_options[lid],loci[lid])
        if cols:write_tsv(prediction_table,raw_predictions,cols)
    write_approved_gff(root,selected,selected_options)
    write_fasta(root/'reviewed_target_cds.fna',[(k,v['cds']) for k,v in sequence_overrides.items()]);write_fasta(root/'reviewed_target_proteins.faa',[(k,v['protein']) for k,v in sequence_overrides.items()])
    write_tsv(root/'reviewed_targets.tsv',selected,['locus_id','group','genome_id','model_id','discovery_model_id','review_option_id','review_origin','target_name','status','window_start','window_end','target_start','target_end','actual_left_flank_bp','actual_right_flank_bp','left_clipped','right_clipped'])
    try:
        finish_comparison(params,source,root,state['original'],state['groups'],selected,state['candidates'],observations,context,proteins,metadata,state['windows'],runner,errors,output_manifest,sequence_overrides)
    except (OSError,ValueError,RuntimeError) as error:errors.append(dict(stage='reviewed_comparison',error=str(error)))
    notice='<aside style="padding:15px;background:#eef2ff"><strong>Manually reviewed inputs.</strong> See <a href="approved_models.tsv">target decisions</a>, <a href="approved_neighbors.tsv">neighbour decisions</a> and <a href="reviewed_targets.tsv">actual target/flank coordinates</a>. Human decisions do not independently prove orthology. Prepared DNA windows are unchanged when the chosen CDS changes.</aside>'
    if (root/'report.html').exists():
        p=root/'report.html';p.write_text(p.read_text().replace('</h1>','</h1>'+notice,1))
    output_manifest.update(finished=now(),commands=runner.commands,errors=errors,status='completed_with_errors' if errors or state['original'].get('errors') else 'completed')
    (root/'manifest.json').write_text(json.dumps(output_manifest,indent=2,default=str)+'\n')
    print('Reviewed comparison: '+str(root),file=sys.stderr)
    return 1 if errors or state['original'].get('errors') else 0


def review_contexts(args):
    """Turn an existing comparative run into STOP 2 without searching/predicting."""
    from .cli import parser
    previous=Path(args.results).resolve();old=json.loads((previous/'manifest.json').read_text())
    source=Path(args.search_results).resolve() if args.search_results else Path(old['source_results'])
    for name,key in [('manifest.json','source_manifest_sha256'),('models.tsv','source_models_sha256')]:
        if old.get(key) and sha256(source/name)!=old[key]:raise ValueError('Imported comparison/search mismatch: '+name)
    original=json.loads((source/'manifest.json').read_text());models={r['model_id']:r for r in table(source/'models.tsv')}
    regions={r['region_id']:r for r in table(source/'dna_regions.tsv')}
    selected=[]
    for record in table(previous/'loci.tsv'):
        model=models.get(record['model_id']) or regions.get(record['model_id'])
        if model is None:raise ValueError('Imported locus absent from source search: '+record['model_id'])
        row={**model,**record}
        for name in ['start','end','target_start','target_end','window_start','window_end','window_bp','softmasked_bases','protein_aa']:
            if row.get(name) not in (None,''):row[name]=int(row[name])
        for name in ['left_clipped','right_clipped']:row[name]=row.get(name) in (True,'True','true')
        if row.get('window_start') is not None:
            row['relative_cds']=[orient_interval(row,a,b) for a,b,*_ in json.loads(row.get('cds_exons') or '[]')]
            row['target_relative_start'],row['target_relative_end']=orient_interval(row,row['target_start'],row['target_end'])
        selected.append(row)
    context=table(previous/'neighbors.tsv')
    for f in context:
        for name in ['start','end','relative_start','relative_end']:f[name]=int(f[name])
        f['clipped']=f['clipped'] in ('True','true')
    metadata={r['protein_id']:{k:v for k,v in r.items() if k!='protein_id'} for r in table(previous/'neighbor_proteins.tsv')}
    protein_path=previous/'neighbor_homology/proteins.faa'
    proteins=optional_fasta(protein_path)
    windows={r['locus_id']:read_fasta(previous/r['window_fasta'],preserve_case=True)[r['locus_id']] for r in selected if not r.get('comparison_error')}
    root=new_root(args.out,previous);_copy_evidence(previous,root)
    config=parser().parse_args(['compare','--results',str(source),'--out',str(root)])
    for key,value in old['parameters'].items():setattr(config,key,value)
    config.results=str(source);config.out=str(root);config._imported_contexts=True
    config._locus_review=dict(status='imported_unreviewed_selection',previous_comparison=str(previous),note='Original automatic loci are not approved. STOP 2 still requires an explicit keep/exclude/unresolved choice for every pair. Alternative loci not present here require STOP 1 and new context extraction.')
    errors=[r for r in old.get('errors',[]) if r.get('stage') in {'window','context','context_inference','context_prediction','context_identification'}]
    manifest={**old,'source_results':str(source),'imported_from':str(previous),'imported_manifest_sha256':sha256(previous/'manifest.json'),'imported_original_errors':old.get('errors',[]),'prediction_run_version':old['version'],'parameters':vars(config)}
    return save_context_review(config,source,root,original,old['groups'],selected,table(previous/'candidates.tsv'),table(previous/'summary.tsv'),context,proteins,metadata,windows,Runner(root),errors,manifest)
