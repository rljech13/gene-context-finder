"""Small derived review kits; native evidence and saved decisions are read-only."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import zipfile

from .sequence import read_fasta
from .alignment_blossom import BLOSSOM

MARKER = '<!-- gcf-interactive-review-20261119-v1 -->'
STYLE = """#interactiveReview{margin:22px 0;padding:20px;border:1px solid #cbd9df;border-radius:8px;background:white;color:#243447;font:14px/1.5 system-ui}#interactiveReview h2{font-size:19px;margin:0 0 8px}.interactive-tools{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:15px 0}.interactive-tools button,.interactive-tools select,.interactive-tools input{font:inherit;padding:7px;border:1px solid #bccbd2;border-radius:5px;background:white}.interactive-tools button{cursor:pointer}.interactive-models{max-height:260px;overflow:auto;display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,360px),1fr));gap:4px 14px;border-top:1px solid #e0e6ea;padding:10px 0}.interactive-model{display:flex;gap:8px;align-items:start;font-size:12px}.interactive-alignment{overflow:auto;max-height:650px;border:1px solid #dce5e8}.interactive-aln-row{display:flex;width:max-content;min-width:100%;font:12px/22px ui-monospace,Consolas,monospace}.interactive-aln-label{position:sticky;left:0;background:#f7fafb;min-width:280px;max-width:280px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;padding:0 8px;z-index:1}.interactive-aln-row code{white-space:pre}.interactive-aln-row code span{display:inline-block;width:9px;text-align:center}#realignStatus{white-space:pre-wrap;font-size:12px}"""

def sequences(root):
    root = Path(root)
    proteins = dict(read_fasta(root/'proteins.faa')) if (root/'proteins.faa').exists() and (root/'proteins.faa').stat().st_size else {}
    coding = dict(read_fasta(root/'cds.fna')) if (root/'cds.fna').exists() and (root/'cds.fna').stat().st_size else {}
    for name, target in [('reconstructed.proteins.faa', proteins), ('reconstructed.cds.fna', coding)]:
        p=root/'reconstruction'/name
        if p.exists(): target.update(read_fasta(p))
    return proteins,coding

def stage2_dataset(root, data):
    protein,cds=sequences(root);records=[];visible=[]
    for m in data['models']:
        mid=m['model_id']
        if mid not in protein or '*' in protein[mid]:continue
        region=next(r for r in data['regions'] if r['region_id']==m['region_id'])
        records.append(dict(id=mid,protein=protein[mid],cds=cds.get(mid,''),group=region.get('gene','target'),
            label=region['genome_id'].split('_GC')[0].replace('_',' ')+' · '+m['original_model']+' · '+m['source']+' · '+str(m['protein_length'])+' aa',
            reference=False,metadata=dict(m,stop1_comment=region.get('notes',''))))
    for region in data['regions']:
        own=[m for m in data['models'] if m['region_id']==region['region_id'] and m['overlaps_anchor'] and m['strand']==region['strand'] and m['model_id'] in protein]
        if own:
            best=min(own,key=lambda m:(m['coding_status']!='complete_predicted',not m.get('target_support'),m['source']!='annotation',m['model_id']))
            visible.append(best['model_id'])
    queries=Path(root)/'target_support/subjects.faa'
    if queries.exists():
        for gene in dict.fromkeys(r.get('gene','target') for r in data['regions']):
            for qid,seq in read_fasta(queries).items():
                rid='reference_'+hashlib.sha256(gene.encode()).hexdigest()[:8]+'_'+qid
                q=(data.get('query_members',{}).get(qid) or [{}])[0]
                records.append(dict(id=rid,protein=seq,cds='',group=gene,reference=True,label=(q.get('species') or 'Query')+' · '+(q.get('protein') or q.get('query') or qid),metadata=dict(query_id=qid,cds_status='not_supplied')))
                visible.append(rid)
    return dict(schema='gcf_interactive_review_v1',stage=2,review_id=data['review_id'],sequences=records,portable_available=False,
        selected=[s['id'] for s in records if s['reference']],visible=visible,colours=BLOSSOM,history=dict(source_review_id=data['review_id'],note='No model approved by opening this package.'))

def stage3_dataset(root, data, source=None, source_display=None):
    """Candidate pool can grow; exact original selected set is kept separately."""
    root=Path(root);raw=data.get('finalization',{}).get('raw_sequences',{})
    if not raw:
        raw={kind:{m['model_id']:data['alignments'][m['group']][kind][m['model_id']].replace('-','') for m in data['models']} for kind in ['protein','cds']}
    records={};selected=[];groups={g['id']:g['name'] for g in data['groups']}
    for m in data['models']:
        mid=m['model_id'];records[mid]=dict(id=mid,protein=raw['protein'][mid],cds=raw['cds'][mid],group=m['group'],
            label=m['species']+' · '+m['original_model']+' · '+m['source'],reference=False,metadata=m);selected.append(mid)
    if source is not None:
        protein,coding=sequences(source)
        source_display=source_display or json.loads((Path(source)/'evidence.json').read_text())
        source_models={m['model_id']:m for m in source_display['models']}
        constraints=source_display.get('model_constraints',{})
        physical=lambda m:constraints.get(m['model_id'],{}).get('physical_id',m['gene_group'])
        group_by_physical={physical(source_models[m['model_id']]):m['group'] for m in data['models'] if m['model_id'] in source_models}
        group_by_region={source_models[m['model_id']]['region_id']:m['group'] for m in data['models'] if m['model_id'] in source_models}
        for m in source_models.values():
            mid=m['model_id']
            gid=group_by_physical.get(physical(m))
            if not gid and m['overlaps_anchor']:gid=group_by_region.get(m['region_id'])
            if not gid or mid in records or mid not in protein or '*' in protein[mid]:continue
            records[mid]=dict(id=mid,protein=protein[mid],cds=coding.get(mid,''),group=gid,reference=False,
                label=m['genome_id'].split('_GC')[0].replace('_',' ')+' · '+m['original_model']+' · '+m['source']+' · alternative CDS',
                metadata=dict(m,comparison_only=True,target_warnings=constraints.get(mid,{}).get('target_warnings',[])))
    for g in data['groups']:
        path=root/'queries'/g['id']/'proteins.input.faa'
        inputs=read_fasta(path) if path.exists() else {}
        for rid,seq in inputs.items():
            if not rid.startswith('reference_'):continue
            members=data.get('query_members',{}).get(rid.removeprefix('reference_'),[]) or [{}]
            q=members[0]
            records[rid]=dict(id=rid,protein=seq,cds='',group=g['id'],reference=True,
                label=(q.get('species') or 'Query')+' · '+(q.get('protein') or q.get('query') or rid),metadata=dict(aliases=members,cds_status='not_supplied'))
            selected.append(rid)
    return dict(schema='gcf_interactive_review_v1',stage=3,review_id=data['review_id'],groups=groups,portable_available=False,
        sequences=list(records.values()),selected=selected,visible=selected,colours=BLOSSOM,
        history=dict(original_stop2_decisions=json.loads((root/'decisions.snapshot.json').read_text()) if (root/'decisions.snapshot.json').exists() else {},
            comments=data.get('comment_history',{}),original_final_choices=data.get('finalization',{}).get('state',{}),
            note='Original choices and reviewed flags are unchanged; interactive membership is not biological approval.'))

def attach(html,data):
    from .synteny_reference import attach as attach_reference
    if MARKER in html:return attach_reference(html)
    text=json.dumps(data,ensure_ascii=True,separators=(',',':')).replace('<',chr(92)+'u003c').replace('&',chr(92)+'u0026')
    section='<section id="interactiveReview"></section><script id="interactiveReviewData" type="application/json">'+text+'</script>'
    if '<details id="alignmentPanel"' in html:
        html=html.replace('<details id="alignmentPanel"',section+'<details id="alignmentPanel"',1)
    else:html=html.replace('</main>',section+'</main>',1)
    html=html.replace('</head>','<script>'+Path(__file__).with_name('review_storage.txt').read_text()+'</script><style>'+STYLE+'</style></head>',1)
    script=Path(__file__).with_name('review_realign.txt').read_text()
    script+='\n'+Path(__file__).with_name('review_undo.txt').read_text()
    from .review_storage import finalize
    return finalize(attach_reference(MARKER+'\n'+html.replace('</body>','<script>'+script+'</script></body>',1)))

def write_dataset(root, data):
    path=Path(root)/'review_data.json'
    path.write_text(json.dumps(data,ensure_ascii=True,separators=(',',':'))+'\n')


def write_kit(root,data,bin_dir=None,platforms=None,*,page=None):
    """The tiny app is independent of source scientific trees and environments."""
    data=dict(data,portable_available=True)
    if page is not None:
        import re
        text=json.dumps(data,ensure_ascii=True,separators=(',',':')).replace('<',chr(92)+'u003c').replace('&',chr(92)+'u0026')
        page,count=re.subn(r'(<script id="interactiveReviewData" type="application/json">).*?(</script>)',lambda m:m[1]+text+m[2],page,count=1,flags=re.S)
        if count!=1:raise ValueError('Expected one embedded review dataset')
    root=Path(root);write_dataset(root,data)
    binary_root=Path(bin_dir or os.environ.get('GCF_REVIEW_BIN_DIR') or Path(__file__).parent/'vendor/review')
    if not binary_root.is_dir() or not list(binary_root.glob('review-*')):
        raise ValueError('Build the optional local launcher with scripts/build_review.py and set GCF_REVIEW_BIN_DIR; ordinary report generation needs no launcher')
    files={'review_data.json':(root/'review_data.json').read_bytes()}
    launcher=r'''#!/bin/sh
cd "$(dirname "$0")" || exit 1
os=$(uname -s | tr '[:upper:]' '[:lower:]')
case "$(uname -m)" in arm64|aarch64) arch=arm64;; x86_64|amd64) arch=amd64;; *) echo "Unsupported CPU architecture"; exit 1;; esac
chmod +x ".review/bin/review-$os-$arch"
exec ".review/bin/review-$os-$arch" --root "$PWD"
'''
    windows='@echo off\r\ncd /d "%~dp0"\r\nset "CPU=amd64"\r\nif /I "%PROCESSOR_ARCHITECTURE%"=="ARM64" set "CPU=arm64"\r\n".review\\bin\\review-windows-%CPU%.exe" --root "%CD%"\r\npause\r\n'
    files.update({'Start_review.command':launcher.encode(),'Start_review.sh':launcher.encode(),'Start_review.cmd':windows.encode()})
    for name,raw in files.items():
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
        if name.endswith(('.command','.sh')):p.chmod(0o755)
    binaries=[p for p in binary_root.glob('review-*') if platforms is None or p.name.removesuffix('.exe').removeprefix('review-') in platforms]
    if not binaries:raise ValueError('No packaged launcher for the requested platform')
    for p in binaries:
        dest=root/'.review/bin'/p.name;dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.is_symlink():dest.unlink()
        if not dest.exists():
            try:os.link(p,dest)
            except OSError:shutil.copyfile(p,dest)
        elif hashlib.sha256(dest.read_bytes()).digest()!=hashlib.sha256(p.read_bytes()).digest():
            raise ValueError('Use a new output directory for changed launcher binaries')
        files['.review/bin/'+p.name]=p.read_bytes()
    mini='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gene context · interactive alignment</title></head><body><main style="max-width:1400px;margin:24px auto;padding:16px"><h1>Gene context · CDS model comparison</h1><p>Exact native models and query references. Original decisions are retained; comparison membership does not approve a target.</p></main></body></html>'
    files['report.html']=(page if page is not None else attach(mini,data)).encode()
    files['README.txt']=b'Extract all files. Start Start_review.command (macOS), Start_review.sh (Linux), or Start_review.cmd (Windows). No Python/Conda/Docker is required. The launcher opens the report on localhost. MAFFT is used from PATH when suitable; otherwise the verified official portable package is downloaded only when Recalculate is pressed. Finish and clean up waits for all calculations and saved results. Other live report sessions retain their leased MAFFT. Native inputs, comments, imported choices and review_results remain. A disconnected tab does not abort a calculation. Unsupported native MAFFT architectures need a preinstalled compatible MAFFT. Source history is not human final approval.\n'
    files['THIRD_PARTY.txt']=b'The local launcher uses the Go standard library. MAFFT is downloaded separately from its official distribution with a pinned SHA256. MAFFT Unix portable packages are BSD; the Windows all-in-one package is GPL and includes MSYS. Upstream sources and notices: https://mafft.ddbj.nig.ac.jp/alignment/software/ . Package notices are retained in the extracted runtime. The residue palette follows Jalview GeCoS Blossom.\n'
    files['LICENSE-Go.txt']=(binary_root/'LICENSE-Go.txt').read_bytes()
    (root/'LICENSE-Go.txt').write_bytes(files['LICENSE-Go.txt'])
    (root/'REVIEW_RUNTIME.txt').write_bytes(files['README.txt']+b'\n'+files['THIRD_PARTY.txt'])
    checks={n:hashlib.sha256(v).hexdigest() for n,v in files.items()}
    files['PACKAGE_CHECKSUMS.json']=json.dumps(checks,indent=2).encode()
    for name,raw in files.items():
        path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists() or path.read_bytes()!=raw:path.write_bytes(raw)
    with zipfile.ZipFile(root/'review_app.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,raw in files.items():
            info=zipfile.ZipInfo('Gene_context_review/'+name);info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(0o100755 if name.startswith('.review/bin/') or name.endswith(('.command','.sh')) else 0o100644)<<16
            z.writestr(info,raw)
    return dict(schema=data['schema'],review_id=data['review_id'],sequences=len(data['sequences']),selected=len(data['selected']),
        app_zip_sha256=hashlib.sha256((root/'review_app.zip').read_bytes()).hexdigest(),app_zip_bytes=(root/'review_app.zip').stat().st_size)


def package_review(args):
    import platform
    from .stage1 import new_output
    source=Path(args.results).resolve();target=Path(args.out).resolve()
    if source==target or source in target.parents or target in source.parents:
        raise ValueError('Use a separate portable-package output directory')
    data=json.loads((source/'review_data.json').read_text())
    if data.get('schema')!='gcf_interactive_review_v1':raise ValueError('Expected saved interactive review data')
    page=None
    if data.get('workflow')=='gcf_alignment_first_v1':
        seal=json.loads((source/'seal.json').read_text())
        for name in ('manifest.json','review_data.json','report.html'):
            if hashlib.sha256((source/name).read_bytes()).hexdigest()!=seal.get(name):
                raise ValueError('Review presentation seal differs: '+name)
        manifest=json.loads((source/'manifest.json').read_text())
        if manifest.get('review_id')!=data['review_id']:
            raise ValueError('Review dataset belongs to another workspace')
        review_id=data['review_id'];page=(source/'report.html').read_text()
    elif data.get('stage')==2:
        from .stage2_choices import verify_checkpoint
        _,scientific=verify_checkpoint(source);review_id=scientific['review_id']
    elif data.get('stage')==3:
        from .stage4 import checked_source
        _,manifest=checked_source(source);review_id=manifest['source_review_id']
    else:raise ValueError('Portable review requires Stage 2 or Stage 3')
    if data['review_id']!=review_id:raise ValueError('Review dataset belongs to another checkpoint')
    chosen=args.platform
    if chosen=='current':
        osname={'Darwin':'darwin','Linux':'linux','Windows':'windows'}.get(platform.system())
        arch={'arm64':'arm64','aarch64':'arm64','x86_64':'amd64','AMD64':'amd64'}.get(platform.machine())
        if not osname or not arch:raise ValueError('Choose a supported viewer platform explicitly')
        chosen=osname+'-'+arch
    root=new_output(target)
    return write_kit(root,data,platforms=None if chosen=='all' else {chosen},page=page)
