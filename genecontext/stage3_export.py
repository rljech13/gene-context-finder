"""Portable default selections with literal STOP-1/2 comments and provenance."""
import csv
import hashlib
from html import escape
import io
import json
from pathlib import Path
import zipfile

from .sequence import sha256, read_fasta


def json_text(value):
    return json.dumps(value, ensure_ascii=False, indent=2)+'\n'


def tsv(rows, fields=None):
    fields=fields or (list(rows[0]) if rows else ['model_id'])
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
    return stream.getvalue()


def fasta(rows, records):
    return ''.join('>'+r['model_id']+'\n'+'\n'.join(records[r['model_id']][i:i+80]
        for i in range(0,len(records[r['model_id']]),80))+'\n' for r in rows)


def collect_history(root, payload):
    """Read saved comments; no reconstruction of absent historical remarks."""
    root=Path(root);evidence=root/'source.evidence.json';digest=sha256(evidence)
    saved=root/'selection_history.json'
    if saved.exists():
        history=json.loads(saved.read_text())
        if history['source_evidence_sha256']!=digest or history['review_id']!=payload['review_id']:
            raise ValueError('Comment history belongs to another STOP-2 checkpoint')
        return history
    data=json.loads(evidence.read_text());state=json.loads((root/'decisions.snapshot.json').read_text())
    regions={r['region_id']:r for r in data['regions']};source_models={m['model_id']:m for m in data['models']}
    comments={};records=[]
    for model in payload['models']:
        mid=model['model_id'];rid=source_models.get(mid,model).get('region_id') or model.get('region_id')
        r=regions[rid]
        rc=state.get('regions',{}).get(rid,{})
        comments[mid]=dict(stop1_reason=r.get('reason',''),stop1_notes=r.get('notes',''),
            stop2_region_notes=rc.get('notes',''),stop2_model_notes=state.get('models',{}).get(mid,{}).get('notes',''),
            stop2_automatic_comments=state.get('automatic_comments',{}).get(rid,[]),
            stop2_review_plan=state.get('review_plan',{}).get(rid,{}))
        for field,value in comments[mid].items():
            values=value if isinstance(value,list) else [value]
            for text in values:
                if text not in ('',{},None):records.append(dict(stage='STOP1' if field.startswith('stop1') else 'STOP2',
                    model_id=mid,genome_id=r['genome_id'],region_id=rid,field=field,
                    text=text if isinstance(text,str) else json.dumps(text,ensure_ascii=False),
                    origin='automatic' if field=='stop2_automatic_comments' else 'saved_decision'))
    # Preserve complete saved decisions too: comments on unchosen alternatives,
    # neighbors and pair choices must not disappear from the portable history.
    files={'comments/stop2.decisions.json':(root/'decisions.snapshot.json').read_bytes().decode('utf-8')}
    provenance=[dict(stage='STOP2',file='comments/stop2.decisions.json',sha256=sha256(root/'decisions.snapshot.json'))]
    manifest=json.loads((root/'manifest.json').read_text()) if (root/'manifest.json').exists() else {}
    locator=manifest.get('source_checkpoint') or payload.get('parameters',{}).get('results')
    source=Path(locator) if locator else None
    status='inherited_region_comments_only'
    if source and (source/'evidence.json').is_file():
        if sha256(source/'evidence.json')!=digest:
            raise ValueError('STOP-2 source differs from the report evidence')
        seal=json.loads((source/'seal.json').read_text()) if (source/'seal.json').exists() else {}
        for name in ['selection.input.json','selection.input.tsv','selection.snapshot.tsv']:
            path=source/name
            if not path.exists():continue
            if name not in seal or sha256(path)!=seal[name]:
                raise ValueError('STOP-1 saved selection is not sealed or changed: '+name)
            if name.startswith('selection.input.') and data.get('selection',{}).get('decision_source_sha256') not in (None,sha256(path)):
                raise ValueError('STOP-1 saved selection fingerprint differs')
            target='comments/stop1.'+name
            files[target]=path.read_bytes().decode('utf-8');provenance.append(dict(stage='STOP1',file=target,sha256=sha256(path)))
            status='saved_input_and_inherited_region_comments'
    return dict(schema='stage3_selection_history_v1',review_id=payload['review_id'],
        source_evidence_sha256=digest,source_stage1_review_id=data.get('source_review_id',''),
        stop1_history_status=status,model_comments=comments,comments=records,files=files,provenance=provenance,
        note='Literal saved comments are retained. Missing historical files are explicit; export does not imply manual approval.')


def model_rows(rows, history):
    return [dict(row,**{k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
        for k,v in history['model_comments'].get(row['model_id'],{}).items()}) for row in rows]


def query_export(root,payload):
    """Export the exact saved query proteins, separately from model decisions."""
    root=Path(root);proteins={};rows=[];seen=set()
    for group in payload.get('groups',[]):
        expected=group.get('queries',[])
        if not expected:continue
        path=root/'queries'/group['id']/'proteins.input.faa'
        if not path.resolve().is_relative_to(root.resolve()):raise ValueError('Unsafe query source path')
        if not path.is_file():raise ValueError('Saved query protein input is missing: '+group['id'])
        raw=read_fasta(path);digest=sha256(path)
        aligned=payload.get('alignments',{}).get(group['id'],{}).get('protein',{})
        for qid in expected:
            rid='reference_'+qid
            if rid not in raw or not raw[rid]:raise ValueError('Saved query protein is missing: '+rid)
            seq=raw[rid]
            if rid in proteins and proteins[rid]!=seq:raise ValueError('Conflicting query protein: '+rid)
            if rid in aligned and aligned[rid].replace('-','')!=seq:raise ValueError('Query alignment differs from saved input: '+rid)
            proteins[rid]=seq
            if (group['id'],rid) in seen:continue
            seen.add((group['id'],rid))
            members=payload.get('query_members',{}).get(qid,[]) or [{}];first=members[0]
            rows.append(dict(reference_id=rid,query_id=qid,group=group['id'],gene=group['name'],
                accession=first.get('protein') or first.get('query') or qid,
                species=first.get('species') or 'Unspecified species',
                symbol=first.get('gene_symbol') or first.get('reference_label') or group['name'],
                origin='Supplied protein query',protein_aa=len(seq),protein_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                cds_status='not_supplied_to_stop3',source_file=str(path.relative_to(root)),source_file_sha256=digest,
                aliases=json.dumps(members,ensure_ascii=False)))
    if set(proteins)&{m['model_id'] for m in payload['models']}:raise ValueError('Query IDs collide with selected model IDs')
    fields=['reference_id','query_id','group','gene','accession','species','symbol','origin','protein_aa','protein_sha256','cds_status','source_file','source_file_sha256','aliases']
    files={'query_references.faa':fasta([dict(model_id=rid) for rid in proteins],proteins),
        'query_references.tsv':tsv(rows,fields),
        'query_references.json':json_text(dict(schema='stage3_query_references_v1',references=rows,
            protein_sequences=len(proteins),cds_status='not_supplied_to_stop3',
            note='Exact saved protein query inputs, with alignment IDs and original accession metadata. Queries are references, not selected genome models. No CDS is inferred from protein.'))}
    return dict(files=files,protein_sequences=len(proteins))


def write_default_export(root,payload,finalization,queries=None):
    root=Path(root);history=collect_history(root,payload)
    (root/'selection_history.json').write_text(json_text(history))
    rows=finalization['rows'];all_rows=rows['primary']+rows['cases'];state=finalization['state']
    ids=[r['model_id'] for r in all_rows]
    if len(set(ids))!=len(ids) or set(ids)!={m['model_id'] for m in payload['models']}:
        raise ValueError('Default export must include every selected STOP-2 model exactly once')
    queries=queries if queries is not None else query_export(root,payload)
    files=dict(history['files'],**queries['files'])
    files.update({'stage3.final.decisions.json':json_text(state),
        'models_with_comments.tsv':tsv(model_rows(all_rows,history)),
        'comments/history.json':json_text({k:v for k,v in history.items() if k!='files'}),
        'comments/comments.tsv':tsv(history['comments'],['stage','model_id','genome_id','region_id','field','text','origin']),
        'README.md':'# Default STOP 3 selection\n\nEvery STOP-2 selected model is retained. primary_targets and special_cases partition that set. models_with_comments.tsv joins the exact model identity and assignment with saved STOP-1/2 comments. comments/ preserves complete available prior decisions, including comments on unchosen alternatives. Quoted TSV cells can contain tabs and newlines; JSON and source snapshots preserve literal text.\n\nThe saved report defaults are exported without per-model clicking. Suggestions remain reviewed=false; downloading does not approve them. Current browser edits belong to the separate Download split package action. No new alignment or model selection is run. Protein/CDS are complete native sequences; terminal CDS stops are retained. stage3.final.decisions.json binds choices to the exact review and sequence-set IDs.\n'})
    for name,part in [('primary_targets',rows['primary']),('special_cases',rows['cases']),('all_selected',all_rows)]:
        files[name+'.tsv']=tsv(part,list(all_rows[0]) if all_rows else None)
        for kind,suffix in [('protein','proteins.faa'),('cds','cds.fna')]:
            raw=finalization['raw_sequences'][kind]
            if set(raw)!=set(ids):raise ValueError('Native sequence set differs from chosen models')
            files[name+'.'+suffix]=fasta(part,raw)
    for name,part in [('primary_targets',rows['primary']),('all_selected',all_rows)]:
        files[name+'_with_queries.proteins.faa']=fasta(part,finalization['raw_sequences']['protein'])+queries['files']['query_references.faa']
    for name in [*queries['files'],'primary_targets_with_queries.proteins.faa','all_selected_with_queries.proteins.faa']:
        (root/name).write_text(files[name],encoding='utf-8')
    files['README.md']+='\nQuery references are included in query_references.faa, with accession, species, aliases and exact saved-input SHA256 in query_references.tsv/json. primary_targets_with_queries.proteins.faa combines the primary representatives with every supplied query; all_selected_with_queries.proteins.faa also retains separate cases. Model-only FASTAs remain available. Query CDS was not supplied to STOP3 and is not inferred.\n'
    content=dict(schema='stage3_default_export_v1',review_id=state['review_id'],set_id=state['set_id'],
        selected_models=len(ids),primary_models=len(rows['primary']),special_case_models=len(rows['cases']),
        query_proteins=queries['protein_sequences'],queries_included=True,
        defaults_only=True,manual_approval_asserted=False,stop1_history_status=history['stop1_history_status'],
        files={name:hashlib.sha256(text.encode()).hexdigest() for name,text in files.items()})
    files['export_manifest.json']=json_text(content)
    package=root/'default_models_with_comments.zip'
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,text in files.items():archive.writestr(name,text.encode())
    record=dict(content,package=package.name,package_sha256=sha256(package))
    (root/'default_export.json').write_text(json_text(record))
    # Only prior comments are added to the browser's current-choice package;
    # its STOP-3 decisions remain the live browser state, never stale defaults.
    return dict(history=history,files=history['files'],record=record)


def write_collection(root,reports):
    """One click for all ready reports, with a separate namespace per gene."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    combined=[];entries=[];files={};keys=set()
    for report in reports:
        key=report['key'];folder=Path(report['results'])
        if key in keys or not key or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in key):
            raise ValueError('Duplicate or unsafe report key')
        keys.add(key);meta=json.loads((folder/'default_export.json').read_text())
        package=folder/meta['package']
        if sha256(package)!=meta['package_sha256']:raise ValueError('Default package changed: '+key)
        with zipfile.ZipFile(package) as archive:
            if archive.testzip():raise ValueError('Corrupt default package: '+key)
            for name,digest in meta['files'].items():
                value=archive.read(name)
                if hashlib.sha256(value).hexdigest()!=digest:raise ValueError('Export member changed: '+name)
            for name in archive.namelist():
                if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError('Unsafe package member')
                files[key+'/'+name]=archive.read(name)
            state=json.loads(archive.read('stage3.final.decisions.json'))
            for row in csv.DictReader(io.StringIO(archive.read('models_with_comments.tsv').decode()),delimiter='\t'):
                combined.append(dict(report_key=key,**row))
        entries.append(dict(key=key,report_url=report['url'],review_id=meta['review_id'],set_id=meta['set_id'],
            selected_models=meta['selected_models'],primary_models=meta['primary_models'],
            special_case_models=meta['special_case_models'],choices=state,package_sha256=meta['package_sha256']))
    fields=list(dict.fromkeys(k for row in combined for k in row)) or ['report_key','model_id']
    files['all_models_with_comments.tsv']=tsv(combined,fields).encode()
    bundle=dict(schema='stage3_default_bundle_v1',defaults_only=True,manual_approval_asserted=False,
        reports=entries,selected_models=len(combined))
    files['stage3.defaults.bundle.json']=json_text(bundle).encode()
    files['README.md']=b'# All ready STOP 3 default selections\n\nOne directory per gene/report. Complete native proteins/CDS, primary and separate-case tables, exact decisions, and all available prior comment snapshots are included. all_models_with_comments.tsv is the combined catalogue. Suggestions are not marked manually reviewed. This package exports the saved defaults; browser edits are exported from each report separately.\n'
    package=root/'all_default_models_with_comments.zip'
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,value in files.items():archive.writestr(name,value)
    (root/'collection.json').write_text(json_text(dict(bundle,package=package.name,package_sha256=sha256(package))))
    rows=''.join('<tr><td><a href="'+escape(e['report_url'],quote=True)+'">'+escape(e['key'])+'</a></td><td>'+str(e['selected_models'])+'</td><td>'+str(e['primary_models'])+'</td><td>'+str(e['special_case_models'])+'</td></tr>' for e in entries)
    (root/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>STOP 3 — export all defaults</title><style>body{font:16px/1.5 system-ui;max-width:1000px;margin:40px auto;padding:24px;color:#213a3d;background:#f3f6f5}a{color:#1d686b}.button{display:inline-block;padding:14px 20px;border-radius:7px;background:#1d686b;color:white;text-decoration:none}td,th{text-align:left;padding:10px;border-bottom:1px solid #dce6e5}table{width:100%;background:white;border-collapse:collapse}</style><h1>STOP 3: all ready genes</h1><p>'+str(len(entries))+' reports · '+str(len(combined))+' selected models</p><p><a class="button" href="'+package.name+'" download>Export all default models + comments · ZIP</a></p><p>One download for all ready genes. Includes primary targets, separate cases, complete native proteins/CDS and saved STOP-1/2/3 comments. No per-gene or per-model saving is needed. Suggestions retain their review status.</p><p>This exports the report defaults. Export browser edits with Download split package inside the relevant report.</p><table><tr><th>Report</th><th>All models</th><th>Primary targets</th><th>Separate cases</th></tr>'+rows+'</table></html>')
    return dict(reports=len(entries),models=len(combined),package_sha256=sha256(package))
