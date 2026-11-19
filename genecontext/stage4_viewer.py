"""Offline 3D presentation of exact structures; never rerun trees or folding."""
import hashlib
import json
from pathlib import Path
import shutil
import time

from .sequence import sha256
from .stage1 import new_output,dump,now
from .stage4_structure import ca_atoms


def vendor_files():
    folder=Path(__file__).with_name('vendor')
    return [folder/name for name in ['3Dmol-min.js','3Dmol-LICENSE.txt','3Dmol-source.json']]


def label_for(data,sid):
    row=next((r for r in data['rows'] if r['model_id']==sid),None)
    if row:return row['species']+' · '+row['original_model'],'Selected model'
    for tree in data['trees']:
        for leaf in tree.get('leaves',[]):
            if sid in leaf['sequence_ids']:
                return leaf['labels'][leaf['sequence_ids'].index(sid)],'Query reference'
    return sid,'Reference'


def write_viewer(root,data):
    root=Path(root); files=vendor_files()
    if not all(p.is_file() for p in files):raise ValueError('Missing pinned local 3Dmol bundle; install the viewer assets first')
    vendor=json.loads(files[2].read_text())
    for path in files[:2]:
        if sha256(path)!=vendor['sha256'][path.name]:raise ValueError('3Dmol asset hash mismatch: '+path.name)
    (root/'vendor').mkdir(exist_ok=True)
    for path in files:shutil.copyfile(path,root/'vendor'/path.name)
    entries=[];by_hash={};aliases={}
    # Prefer native selected models as representatives over identical query aliases.
    selected={r['model_id'] for r in data['rows']}
    for sid,s in sorted(data['structures'].items(),key=lambda item:item[0] not in selected):
        key=s['protein_sha256'];label,kind=label_for(data,sid)
        path=root/(s.get('full_file') or s['file'])
        if sha256(path)!=s.get('full_file_sha256',s['file_sha256']):raise ValueError('Viewer coordinates changed: '+sid)
        atoms=ca_atoms(path)
        if hashlib.sha256(''.join(a['aa'] for a in atoms).encode()).hexdigest()!=key:raise ValueError('Viewer sequence differs: '+sid)
        entry=by_hash.get(key)
        if not entry:
            entry=dict(id=sid,label=label,kind=kind,source=s['source'],accession=s.get('accession',''),
                file=str(path.relative_to(root)),pdb=path.read_text(),length=s['full_length'],
                coordinate_scope=s.get('coordinate_scope','CA_trace'),
                confidence=s.get('confidence_kind','not_available'),aliases=[])
            entries.append(entry);by_hash[key]=entry
        entry['aliases'].append(dict(id=sid,label=label,kind=kind))
        aliases[sid]=entry['id'];s['viewer_file']='structure_viewer.html'
    payload=dict(entries=entries,aliases=aliases,library=dict(name='3Dmol.js',version=vendor['version']),
        comparisons=data.get('structure_similarity',[]),
        defaults={r['model_id']:r['structure_match']['reference_id'] for r in data['rows'] if r.get('structure_match')},
        note='Exact-sequence structures. Display is not functional or orthology validation.')
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    (root/'structure_viewer.html').write_text(Path(__file__).with_name('stage4_viewer.txt').read_text().replace('__VIEWER_DATA__',text))
    return dict(unique_structures=len(entries),selected_models=sum(sid in selected for sid in aliases),
        viewer_file='structure_viewer.html',library_version=vendor['version'],library_sha256=vendor['sha256']['3Dmol-min.js'])


def attach(args):
    started=time.perf_counter();source=Path(args.results).resolve();target=Path(args.out).resolve()
    if source==target or source in target.parents or target in source.parents:raise ValueError('Use a separate viewer presentation directory')
    seal=json.loads((source/'seal.json').read_text())
    for name,digest in seal.items():
        path=(source/name).resolve()
        if source not in path.parents or not path.is_file() or sha256(path)!=digest:raise ValueError('STOP 4 seal mismatch: '+name)
    old_manifest=json.loads((source/'manifest.json').read_text())
    if old_manifest.get('stage') not in ['stage4_prototype','stage4_viewer_presentation'] or old_manifest['status']!='completed':raise ValueError('Completed STOP 4 required')
    data=json.loads((source/'evidence.json').read_text())
    catalogue=json.loads(Path(args.structures).read_text()) if args.structures else None
    if catalogue and catalogue.get('schema')!='stage4_structures_v1':raise ValueError('Expected stage4_structures_v1')
    replacements={}
    for entry in catalogue['entries'] if catalogue else []:
        sid=entry['sequence_id']
        if sid not in data['structures'] or sid in replacements:raise ValueError('Unknown or duplicate viewer structure: '+sid)
        old=data['structures'][sid]
        path=(Path(args.structures).resolve().parent/entry['file']).resolve()
        if sha256(path)!=entry['file_sha256'] or entry['protein_sha256']!=old['protein_sha256']:raise ValueError('Full-atom structure checksum differs: '+sid)
        full_ca=ca_atoms(path,entry.get('chain','A'));saved_ca=ca_atoms(source/old['file'])
        def signature(atoms):return [(a['aa'],a['line'][30:54],a['confidence']) for a in atoms]
        if signature(full_ca)!=signature(saved_ca):raise ValueError('Full-atom CA coordinates differ from the measured structure: '+sid)
        replacements[sid]=(path,entry)
    root=new_output(target)
    for name in seal:
        dest=root/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,dest)
    for sid,(path,entry) in replacements.items():
        s=data['structures'][sid];dest=root/'structures'/(sid+'.full.pdb')
        atoms=[line for line in path.read_text().splitlines() if line.startswith('ATOM  ') and line[21:22]==entry.get('chain','A') and line[16:17] in (' ','A')]
        dest.write_text('\n'.join(atoms)+'\nEND\n')
        s.update(full_file=str(dest.relative_to(root)),full_file_sha256=sha256(dest),
                 coordinate_scope='all_atom_chain' if len(atoms)>s['full_length'] else 'CA_trace')
    for r in data['rows']:
        r['structure']=data['structures'].get(r['model_id'])
    geometry=None
    if getattr(args,'usalign',None):
        from .stage4_similarity import build_similarity
        geometry=build_similarity(root,data,args.usalign)
    data['structure_viewer']=write_viewer(root,data)
    dump(root/'evidence.json',data)
    text=json.dumps(data,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    (root/'report.html').write_text(Path(__file__).with_name('stage4_report.txt').read_text().replace('__DATA__',text))
    shutil.copyfile(source/'manifest.json',root/'source.stage4.manifest.json')
    dump(root/'manifest.json',dict(stage='stage4_viewer_presentation',status='completed',created=now(),
        source_manifest_sha256=sha256(source/'manifest.json'),source_seal_sha256=sha256(source/'seal.json'),
        structures_catalogue_sha256=sha256(args.structures) if args.structures else None,
        source_code_sha256={n:sha256(Path(__file__).with_name(n)) for n in ['stage4_viewer.py','stage4_viewer.txt','stage4_report.txt','stage4_similarity.py']},
        viewer=data['structure_viewer'],elapsed_seconds=time.perf_counter()-started,
        scientific_recalculation=bool(geometry),structural_geometry=geometry,tree_or_sequence_recalculation=False,inherited_structure_similarity=bool(data.get('structure_similarity') and not geometry),biological_approval_asserted=False,preview=data['preview']))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    print(json.dumps(dict(report=str(root/'report.html'),viewer=data['structure_viewer'],elapsed_seconds=time.perf_counter()-started)))
    return 0
