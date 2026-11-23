"""Rebuild completed reports in the current review form, without new inference.

The original result files and review IDs remain evidence. A deterministic new
review namespace isolates editable drafts from the original completed study.
"""
import argparse
import base64
from copy import deepcopy
import gzip
import hashlib
import html as html_module
import json
from pathlib import Path
import re
import uuid
import zipfile
import shutil
import os

ASSETS = Path(__file__).with_name('vendor')/'legacy_review'


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return json.dumps(value, ensure_ascii=True, separators=(',', ':')).replace('<', chr(92)+'u003c').replace('&', chr(92)+'u0026')


def read_payload(path):
    text = Path(path).read_text()
    if 'const DATA=' not in text:
        for match in re.finditer(r"atob\('([A-Za-z0-9+/=]+)'\)", text):
            try:
                candidate = gzip.decompress(base64.b64decode(match[1])).decode()
            except Exception:
                continue
            if 'const DATA=' in candidate:
                text = candidate
                break
    offset = text.index('const DATA=')+len('const DATA=')
    value, _ = json.JSONDecoder().raw_decode(text[offset:])
    return value


def backtranslate(row, member):
    native = member['cds']
    if len(native) == (len(member['protein'])+1)*3:
        native = native[:-3]
    if len(native) != len(member['protein'])*3 or row.replace('-', '') != member['protein']:
        raise ValueError('Native protein/CDS length or saved alignment differs: '+member['id'])
    codons = iter(native[i:i+3] for i in range(0, len(native), 3))
    return ''.join('---' if aa == '-' else next(codons) for aa in row)


def import_analysis(group, alignment, members, review_id, manifest):
    protein = deepcopy(alignment['protein'])
    if len(protein) < 2 or len({len(row) for row in protein.values()}) != 1:
        raise ValueError('Invalid saved protein alignment: '+group)
    for mid, row in protein.items():
        if mid not in members or row.replace('-', '') != members[mid]['protein']:
            raise ValueError('Saved alignment does not match the exact native member: '+mid)
    cds = {mid:backtranslate(row, members[mid]) for mid, row in protein.items() if members[mid].get('cds')}
    identity = sha(encoded([review_id, group, protein, cds]).encode())
    common = dict(schema='gcf_review_result_v1', review_id=review_id, group=group,
        ids=list(protein), protein=protein, cds=cds, status='done',
        mafft_version=manifest.get('tools', {}).get('mafft', {}).get('version', 'Saved original alignment'),
        provenance=dict(legacy_import=True, new_alignment_calculation=False,
                        source_manifest_sha256=sha(encoded(manifest).encode())))
    msa = dict(common, id=identity[:20], phase='alignment')
    aliases = {}
    for mid in protein:
        digest = sha(members[mid]['protein'].encode())
        aliases.setdefault(digest, []).append(mid)
    leaves = [dict(tree_id='protein_'+digest[:16], protein_sha256=digest,
                   sequence_ids=ids, included=True) for digest, ids in aliases.items()]
    # This is an unavailable-result record, never an invented topology. Original
    # STOP4 trees, where present, remain linked as historical evidence.
    tree = dict(common, id=sha((identity+'unavailable-original-tree').encode())[:20],
        phase='tree', alignment_id=msa['id'],
        phylogeny=dict(status='not_computed',
            reason='No tree bound to this exact saved alignment was imported. No new tree was calculated during the report rebuild.',
            leaves=leaves, outgroup_ids=[], root_basis='', unrooted_newick=None,
            rooted_newick=None, unrooted_svg=None, rooted_svg=None))
    return dict(alignment=msa, tree=tree)


def render(payload, pool, migration, title):
    assets = {p.name:p.read_text() for p in ASSETS.glob('*.txt')}
    template = assets['stage2_review.txt']
    for marker, name in [('__NEIGHBORHOOD__','stage2_neighborhood.txt'), ('__CLUSTER__','stage2_cluster.txt'),
                         ('__REVIEW_QUEUE__','stage2_review_queue.txt'), ('__MODEL_COMPARE__','stage2_model_compare.txt')]:
        template = template.replace(marker, assets[name])
    page = template.replace('__DATA__', encoded(payload))
    page = page.replace('</head>', '<meta name="gcf-review-flow" content="sequential-v2"><meta name="gcf-portable-review" content="review23-v2"><style>'+assets['interactive_style.txt']+'</style></head>', 1)
    controls = '<section id="interactiveReview"></section><script id="interactiveReviewData" type="application/json">'+encoded(pool)+'</script>'
    page = page.replace('</main>', controls+'</main>', 1)
    banner = '<section id="legacyStudyNotice"><h2>'+html_module.escape(title)+' · saved analysis in the current review form</h2><p>Original models, comments and alignment results are retained. Imported proposals keep their original review status. No new alignment, tree, annotation or folding calculation was started.</p><p><a href="original_report.html">Original completed report</a> · <a href="presentation_rebuild.json">Import provenance</a> · <a href="review_app.zip">Download complete local review kit</a></p><p id="legacyImportNotice"></p><label>Import earlier model/final choices <input id="legacyChoiceImport" type="file" accept=".json"></label></section>'
    page = page.replace('<main>', '<main>'+banner, 1)
    bootstrap = '<script id="legacyReviewData" type="application/json">'+encoded(migration)+'</script><script>'+ASSETS.joinpath('legacy_bootstrap.txt').read_text()+'</script>'
    page = page.replace("<script>\n'use strict';const DATA=", bootstrap+"<script>\n'use strict';const DATA=", 1)
    if bootstrap not in page:
        raise ValueError('Native review script boundary changed')
    scripts = [assets['review_realign.txt'], assets['alignment_review_steps.txt'], assets['review_undo.txt'],
               assets['synteny_review_page.txt'], assets['stop4_handoff.txt'], ASSETS.joinpath('legacy_controls.txt').read_text(), assets['synteny_reference_approval.txt']]
    scripts[1] = scripts[1].replace("'DNMT3B · choose models → alignment → tree'", json.dumps(title)+"+' · choose models → alignment → tree'")
    scripts[1] = scripts[1].replace("'DNMT3B_'+", json.dumps(title+'_')+'+')
    # Missing historical trees remain explicitly unavailable, even though their
    # integrity/provenance record can be restored alongside the saved alignment.
    scripts[1] = scripts[1].replace("'Tree verified and saved · '+r.id", "(r.phylogeny?.status==='computed'?'Tree verified and saved · ':'Tree not previously available · ')+r.id")
    page = page.replace('</body>', ''.join('<script>'+s+'</script>' for s in scripts)+'</body>', 1)
    css = ASSETS.joinpath('link_contrast.css').read_text()
    page = page.replace('</head>', '<style id="gcf-synteny-link-contrast">'+css+'</style></head>', 1)
    from .review_storage import finalize
    return finalize(page)


def build(source, out, title=None, confirmed_decisions=None):
    source, out = Path(source), Path(out)
    if out.exists():
        raise ValueError('Use a new presentation directory; original reports are preserved')
    manifest = json.loads((source/'manifest.json').read_text())
    prepared = Path(manifest['source_checkpoint'])
    native_path = prepared/'report.html'
    native_seal = json.loads((prepared/'seal.json').read_text())
    if native_seal.get('report.html') != sha(native_path.read_bytes()):
        raise ValueError('Original prepared model report differs from its seal')
    original = read_payload(native_path)
    legacy = read_payload(source/'report.html')
    consumed = json.loads((source/'decisions.snapshot.json').read_text())
    old_pool = json.loads((source/'review_data.json').read_text())
    if consumed['review_id'] != original['review_id'] or legacy['review_id'] != original['review_id']:
        raise ValueError('Completed comparison and prepared model decisions do not belong together')
    old_state = deepcopy(legacy['finalization']['state'])
    if confirmed_decisions:
        confirmed = json.loads(Path(confirmed_decisions).read_text())
        if (confirmed.get('review_id'), confirmed.get('set_id')) != (old_state['review_id'], old_state['set_id']):
            raise ValueError('Confirmed final choices belong to another original sequence set')
        old_state = confirmed
    title = title or ', '.join(g['name'] for g in legacy['groups'])
    names = {g['id']:g['name'] for g in legacy['groups']}
    source_id = str(uuid.uuid5(uuid.NAMESPACE_URL, 'gcf-legacy-modern-v1|'+original['review_id']+'|'+old_state['set_id']))
    review_id = str(uuid.uuid5(uuid.NAMESPACE_URL, source_id+'|alignment-review'))
    payload = deepcopy(original)
    payload['review_id'] = source_id
    payload['source_review_id'] = original['review_id']
    if payload.get('annotation_context_external'):
        supplement = prepared/'annotation_context.json'
        payload['annotation_features'] = json.loads(supplement.read_text())['features']
        payload['annotation_context_external'] = False
    native_models = {m['model_id']:m for m in payload['models']}
    sequences = []
    for entry in old_pool['sequences']:
        s = deepcopy(entry)
        s['group'] = names.get(s['group'], s['group'])
        if not s['reference']:
            native = native_models[s['id']]
            if sha(s['protein'].encode()) != native['protein_sha256'] or sha(s['cds'].encode()) != native['cds_sha256']:
                raise ValueError('Exact native protein/CDS differs: '+s['id'])
            s['metadata'] = dict(s.get('metadata', {}), **deepcopy(native))
            s['metadata']['model_constraints'] = deepcopy(payload.get('model_constraints', {}).get(s['id'],
                dict(physical_id=native['gene_group'], target_warnings=(["Outside selected interval"] if not native.get('overlaps_anchor') else [])+(["No retained target-query support"] if not native.get('target_support') else []))))
        sequences.append(s)
    selected = {m['model_id'] for m in legacy['models']}
    consumed_targets = {mid for mid,c in consumed['models'].items() if c['decision']=='inspect' and c['role'] in ('target','alignment_candidate')}
    if selected != consumed_targets:
        raise ValueError('Consumed model set differs from completed comparison; explicit mapping is required')
    baseline = deepcopy(consumed)
    baseline['review_id'] = source_id
    pool = dict(old_pool, schema='gcf_interactive_review_v1', workflow='gcf_alignment_first_v1', stage=2,
        review_id=review_id, source_review_id=source_id, sequences=sequences,
        selected=[s['id'] for s in sequences if s['id'] in selected or s['reference']],
        visible=[s['id'] for s in sequences if s['id'] in selected or s['reference']],
        legacy_import=True, history=dict(kind='completed_legacy_rebuild', original_review_id=original['review_id'],
            original_set_id=old_state['set_id'], no_new_scientific_approval=True))
    members = {s['id']:s for s in sequences}
    analyses = {names[gid]:import_analysis(names[gid], aln, members, review_id, manifest) for gid,aln in legacy['alignments'].items()}
    migration = dict(schema='gcf_legacy_review_import_v1', original_review_id=original['review_id'],
        source_review_id=source_id, alignment_review_id=review_id, original_set_id=old_state['set_id'],
        baseline=baseline, analyses=analyses, default_final=old_state,
        provenance=dict(original_review_id=original['review_id'],original_set_id=old_state['set_id'],consumed_decisions_sha256=sha((source/'decisions.snapshot.json').read_bytes()),new_scientific_approval=False),
        buckets=legacy['finalization']['spec']['buckets'], members=[s for s in sequences if s['id'] in selected],
        first_region=native_models[legacy['models'][0]['model_id']]['region_id'])
    page = render(payload, pool, migration, title)
    out.mkdir(parents=True)
    (out/'report.html').write_text(page)
    (out/'review_data.json').write_text(json.dumps(pool,ensure_ascii=True,separators=(',',':'))+'\n')
    (out/'context.html').write_text('<!doctype html><meta charset="utf-8"><title>STOP 3 · synteny</title><p><a href="report.html?view=synteny">Open STOP 3 · synteny and final models</a></p><script>const u=new URL("report.html",location.href);u.searchParams.set("view","synteny");for(const[k,v]of new URLSearchParams(location.search))u.searchParams.set(k==="model"?"legacy_model":k,v);location.replace(u.href);</script>')
    for part in analyses.values():
        for result in part.values():
            folder=out/'review_results'/result['id'];folder.mkdir(parents=True)
            raw=(json.dumps(result,ensure_ascii=False,separators=(',',':'))+'\n').encode();(folder/'result.json').write_bytes(raw)
            (folder/'checksums.json').write_text(json.dumps({'result.json':sha(raw)})+'\n')
    receipt=dict(schema='gcf_completed_review_rebuild_v1',original_review_id=original['review_id'],
        original_set_id=old_state['set_id'],derived_source_review_id=source_id,alignment_review_id=review_id,
        source_manifest_sha256=sha((source/'manifest.json').read_bytes()),source_report_sha256=sha((source/'report.html').read_bytes()),
        consumed_decisions_sha256=sha((source/'decisions.snapshot.json').read_bytes()),
        confirmed_decisions_sha256=sha(Path(confirmed_decisions).read_bytes()) if confirmed_decisions else None,
        original_models=len(selected),query_references=sum(s['reference'] for s in sequences),
        imported_alignments=len(analyses),new_alignment_calculations=0,new_tree_calculations=0,
        original_model_sequence_hashes={mid:{k:sha(members[mid][k].encode()) for k in ('protein','cds')} for mid in sorted(selected)},
        original_reviewed_assignments=sum(c.get('reviewed') is True for c in old_state['choices'].values()),
        new_human_approval_asserted=False,original_files_changed=False)
    (out/'presentation_rebuild.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def package(source, out, runtime, archive, historical=None):
    """Overlay a rebuilt page on unchanged scientific files, then seal the kit.

    Source files and prior archives are never rewritten. All archive members are
    verified after writing; the local app includes the cached analysis results.
    """
    source, out, runtime, archive = map(Path, (source, out, runtime, archive))
    for name in ('report.html', 'context.html', 'review_data.json', 'seal.json',
                 'portable.seal.json', 'portable.checksums.json', 'presentation.json'):
        if (source/name).exists() and not (out/('original_'+name)).exists():
            (out/('original_'+name)).symlink_to((source/name).resolve())
    for p in source.iterdir():
        if not (out/p.name).exists() and p.name not in ('seal.json', 'portable.seal.json', 'portable.checksums.json', '.review', 'review_app.zip'):
            (out/p.name).symlink_to(p.resolve())
    # Replace presentation launchers only, retaining their old bytes separately.
    launcher = '''#!/bin/sh
cd "$(dirname "$0")" || exit 1
os=$(uname -s | tr '[:upper:]' '[:lower:]')
case "$(uname -m)" in arm64|aarch64) arch=arm64;; x86_64|amd64) arch=amd64;; *) echo "Unsupported CPU architecture"; exit 1;; esac
chmod +x ".review/bin/review-$os-$arch"
exec ".review/bin/review-$os-$arch" --root "$PWD"
'''
    def write(name, raw):
        p=out/name
        if p.is_symlink(): p.unlink()
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
    for name in ('Start_review.command','Start_review.sh'):
        write(name,launcher.encode());(out/name).chmod(0o755)
    write('Start_review.cmd',b'@echo off\r\ncd /d "%~dp0"\r\nset "CPU=amd64"\r\nif /I "%PROCESSOR_ARCHITECTURE%"=="ARM64" set "CPU=arm64"\r\n".review\\bin\\review-windows-%CPU%.exe" --root "%CD%"\r\npause\r\n')
    for p in (runtime/'runtimes').glob('review-*'):
        dest=out/'.review/bin'/p.name;dest.parent.mkdir(parents=True,exist_ok=True);
        if not dest.exists():dest.symlink_to(p.resolve())
    wheel=next(runtime.glob('*.whl'))
    if not (out/wheel.name).exists():(out/wheel.name).symlink_to(wheel.resolve())
    write('Build_STOP4.py',(runtime/'Build_STOP4.py').read_bytes())
    write('START_HERE.html',b'<!doctype html><html lang="en"><meta charset="utf-8"><title>Saved study in the current review form</title><h1>Gene context review</h1><p><a href="report.html">Open the complete report</a> &middot; <a href="report.html?view=synteny">STOP 3: synteny and final models</a></p><p>For editing and recalculation, extract the archive and run Start_review for your operating system. Existing decisions and cached results are retained. Use Build_STOP4.py only after reviewing and exporting the STOP 4 input.</p></html>')
    if historical and Path(historical).is_dir():
        directory=out/'historical_trees';directory.mkdir(exist_ok=True)
        links=[]
        for p in sorted(Path(historical).rglob('*')):
            if p.is_file():
                rel=p.relative_to(historical);dest=directory/rel;dest.parent.mkdir(parents=True,exist_ok=True);
                if not dest.exists():dest.symlink_to(p.resolve())
                if p.name=='tree.svg':links.append('<p><a href="'+html_module.escape(str(rel))+'">'+html_module.escape(str(rel.parent))+'</a></p><object type="image/svg+xml" data="'+html_module.escape(str(rel))+'" style="width:100%;height:850px"></object>')
        write('historical_trees/index.html',('<!doctype html><html lang="en"><meta charset="utf-8"><h1>Previously computed STOP 4 trees</h1><p>Preserved original results. Their original inputs and evidence are included here. These trees were not recalculated or reassigned to the imported STOP 2 alignment.</p>'+''.join(links)+'</html>').encode())
        page=(out/'report.html').read_text().replace('<p id="legacyImportNotice">','<p><a href="historical_trees/index.html">View preserved trees from the completed study</a></p><p id="legacyImportNotice">')
        if 'View preserved trees from the completed study' not in (out/'report.html').read_text():write('report.html',page.encode())
    write('presentation.json',(json.dumps({'schema':'gcf_completed_review_presentation_v1','form':'sequential-v2','legacy_import':True,'scientific_recomputation':False,'manual_choices_changed':False,'provenance':'presentation_rebuild.json'},indent=2)+'\n').encode())
    # A compact standalone local kit and a full report archive share one UI.
    files={}
    for directory, _, names in os.walk(out,followlinks=True):
        for name in names:
            p=Path(directory)/name;rel=str(p.relative_to(out))
            if rel in ('review_app.zip','PACKAGE_CHECKSUMS.json','seal.json','portable.seal.json','portable.checksums.json'):continue
            if p.is_file():files[rel]=p
    checks={name:sha(p.read_bytes()) for name,p in files.items()}
    write('PACKAGE_CHECKSUMS.json',(json.dumps(checks,indent=2)+'\n').encode());files['PACKAGE_CHECKSUMS.json']=out/'PACKAGE_CHECKSUMS.json'
    def zipped(path, items, prefix=''):
        with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as z:
            for name,p in sorted(items.items()):z.write(p,prefix+name)
        with zipfile.ZipFile(path) as z:
            if len(z.namelist())!=len(set(z.namelist())):raise ValueError('Duplicate archive members')
            for name,p in items.items():
                if sha(z.read(prefix+name))!=sha(p.read_bytes()):raise ValueError('Archive member differs: '+name)
    zipped(out/'review_app.zip',files)
    files['review_app.zip']=out/'review_app.zip'
    seal={name:sha(p.read_bytes()) for name,p in files.items()}
    write('seal.json',(json.dumps(seal,indent=2)+'\n').encode());files['seal.json']=out/'seal.json'
    archive.parent.mkdir(parents=True,exist_ok=True);zipped(archive,files,source.name+'/')
    unchanged={str(p.relative_to(source)):sha(p.read_bytes()) for directory,_,names in os.walk(source,followlinks=True) for name in names if (p:=Path(directory)/name).is_file() and str(p.relative_to(source)) in checks and checks[str(p.relative_to(source))]==sha(p.read_bytes())}
    return dict(archive_sha256=sha(archive.read_bytes()),archive_bytes=archive.stat().st_size,
                verified_members=len(files),unchanged_original_files=len(unchanged),unchanged_sha256=unchanged)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--title');p.add_argument('--confirmed-decisions',type=Path);a=p.parse_args()
    print(json.dumps(build(a.results,a.out,a.title,a.confirmed_decisions)))


if __name__=='__main__':
    main()
