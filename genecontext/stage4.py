"""STOP 4 prototype: a compact, traceable candidate dossier from sealed STOP 3."""
from collections import Counter
import hashlib
import importlib
import json
from pathlib import Path
import re
import shutil
import time

from .sequence import read_fasta, write_fasta, write_tsv, sha256, translate
from .stage1 import new_output, dump, now
from .stage3_analysis import read_table
from .stage3_finalization import validate_choices, rows_for
from .stage3 import TimedRunner
from .stage4_phylogeny import build_tree
from .stage4_structure import import_structures, compare_structures


def checked_source(folder):
    root = Path(folder).resolve()
    seal = json.loads((root/'seal.json').read_text())
    required = {'manifest.json','comparison.json','sequence_catalogue.tsv','selected.proteins.faa',
                'selected.cds.fna','finalization.evidence.json','finalization.decisions.json',
                'reference_comparisons.tsv','model_quality.tsv','synteny.tsv'}
    if not required <= set(seal): raise ValueError('STOP 3 seal lacks required evidence files')
    for name,digest in seal.items():
        path = (root/name).resolve()
        if root not in path.parents or not path.is_file() or sha256(path) != digest:
            raise ValueError('STOP 3 seal mismatch: '+name)
    manifest = json.loads((root/'manifest.json').read_text())
    if manifest.get('stage') != 'stage3_comparison' or manifest.get('status') != 'completed' or manifest.get('errors'):
        raise ValueError('A completed, error-free STOP 3 checkpoint is required')
    return root, manifest


def load_models(source, comparison, evidence):
    catalogue = read_table(source/'sequence_catalogue.tsv')
    quality = {q['model_id']:q for q in read_table(source/'model_quality.tsv')}
    rows = {r['model_id']:r for r in catalogue}
    selected = {r['model_id']:r for r in comparison['selected']}
    proteins = read_fasta(source/'selected.proteins.faa')
    coding = read_fasta(source/'selected.cds.fna')
    if len(rows) != len(catalogue) or set(rows) != set(selected) or set(rows) != set(proteins) or set(rows) != set(coding) or set(rows) != set(evidence['models']):
        raise ValueError('Selected model IDs differ across STOP 3 inputs')
    groups = {g['id']:g['name'] for g in comparison['groups']}
    models = []
    if any(not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]{0,150}', x) for x in [*rows,*groups]):
        raise ValueError('Unsafe exported model or group ID')
    for mid,r in rows.items():
        for field,seq in [('protein_sha256',proteins[mid]),('cds_sha256',coding[mid])]:
            if hashlib.sha256(seq.encode()).hexdigest() != r[field]: raise ValueError('Sequence hash mismatch: '+mid)
        if translate(coding[mid]).removesuffix('*') != proteins[mid] or '*' in proteins[mid]:
            raise ValueError('CDS/protein translation differs: '+mid)
        ints = {k:int(r[k]) for k in ['start','end','protein_aa']}
        if ints['protein_aa'] != len(proteins[mid]) or ints['start'] >= ints['end']:
            raise ValueError('Model length or genomic interval differs: '+mid)
        models.append(dict(r,**ints,group=selected[mid]['group'],quality=quality[mid],
            target_gene=groups[selected[mid]['group']],**evidence['models'][mid]))
    signature = [(m['group'],m['genome_id'],m['model_id'],m['cds_sha256']) for m in models]
    expected = hashlib.sha256(json.dumps([evidence['spec']['review_id'],sorted(signature)],sort_keys=True).encode()).hexdigest()
    if expected != evidence['spec']['set_id']: raise ValueError('Finalization selected-set signature differs')
    buckets = evidence['spec']['buckets']
    ids = [mid for b in buckets for mid in b['model_ids']]
    if len(ids) != len(set(ids)) or set(ids) != set(rows): raise ValueError('Final buckets do not partition selected models')
    for b in buckets:
        if any(selected[mid]['group'] != b['group'] or selected[mid]['genome_id'] != b['genome_id'] for mid in b['model_ids']):
            raise ValueError('Foreign model in final gene/genome bucket')
    return models,proteins,coding


def assessment(model, reference_rows, synteny):
    mid = model['model_id']; assignment = model['assignment']; q = model['quality']
    best = assignment.get('best_target') or {}
    hit = next((r for r in reference_rows if r['model_id']==mid and r['kind']=='protein' and r['reference_id']==best.get('reference_id')),None)
    facts,limitations = [], []
    if hit:
        facts.append(f"Target reference {hit['reference_accession']}: identity {float(hit['identity']):.1%}; model/reference coverage {float(hit['coverage_model']):.1%} / {float(hit['coverage_reference']):.1%}.")
    else: limitations.append('No measured comparison to a target query.')
    ident = model['copy_identity']
    facts.append('Source copy label: '+ident['display_name']+' ('+ident['basis'].replace('_',' ')+').')
    conflict = assignment.get('known_other_copy') or assignment.get('alternate_preferred')
    if conflict: identity = 'Conflicting copy evidence'
    elif hit and min(float(hit['coverage_model']),float(hit['coverage_reference'])) >= .5 and ident['status'] not in {'unresolved','homology_only'}:
        identity = 'Consistent candidate'
    else: identity = 'Identity unresolved'
    limitations.append('Specific orthology is not established: no automatic target-versus-paralog clade assignment was made.')
    anchors = [s for s in synteny if mid in [s['locus_a'],s['locus_b']] and int(s.get('unique_anchors') or 0)>=2]
    if anchors: facts.append(f'{len(anchors)} comparisons have at least two reviewed neighbor anchors; this is context evidence, not a gene-name assignment.')
    else: limitations.append('Synteny unresolved: insufficient reviewed neighbor anchors.')
    problems = []
    if q['coding_status'] != 'complete_predicted': problems.append('Partial coding model')
    if str(model.get('window_edge','')).lower() == 'true': problems.append('Saved window edge')
    if q['source_agreement']=='different_CDS_models': problems.append('Annotation / Tiberius differ')
    if q.get('cautions'): limitations.append('STOP 3 flags: '+q['cautions'])
    facts.append('Saved protein/CDS hashes and CDS translation verified.')
    limitations.append('Transcript existence and expression are not evaluated; prediction or an accession alone is not experimental mRNA validation.')
    return dict(identity_assessment=identity,model_assessment='Check model' if problems else 'Coding consistent',
        model_notes='; '.join(problems) or 'Native coding-sequence checks passed',
        synteny_status='Reviewed anchors available' if anchors else 'Not resolved',
        reviewed_anchor_comparisons=len(anchors),facts=facts,limitations=limitations,
        reference_accession=hit['reference_accession'] if hit else '',
        protein_identity=float(hit['identity']) if hit else None,
        coverage_model=float(hit['coverage_model']) if hit else None,
        coverage_reference=float(hit['coverage_reference']) if hit else None)


def reference_panel(path, groups):
    if not path: return []
    data = json.loads(Path(path).read_text())
    if data.get('schema') != 'stage4_reference_panel_v1': raise ValueError('Expected stage4_reference_panel_v1')
    by_name = {g['name']:g['id'] for g in groups}; allowed={g['id'] for g in groups}
    entries=[]; seen=set()
    for r in data['entries']:
        r=dict(r); gid=by_name.get(r['group'],r['group']); sid=r['id']
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,79}',sid) or sid in seen: raise ValueError('Invalid or duplicate extra-reference ID')
        seen.add(sid)
        if gid not in allowed or r['role'] not in ['target','paralog','outgroup'] or not r.get('source'):
            raise ValueError('Extra reference needs a known group, role and source')
        if not r['sequence'] or set(r['sequence'])-set('ACDEFGHIKLMNPQRSTVWY'): raise ValueError('Extra reference must contain an unaligned canonical protein sequence')
        if r['role']=='outgroup' and not r.get('root_basis','').strip(): raise ValueError('Outgroup requires a biological root_basis')
        r['group']=gid;entries.append(r)
    return entries


def run(args):
    if args.bootstrap<1000: raise ValueError('UFBoot requires at least 1000 replicates')
    started=time.perf_counter(); source,manifest=checked_source(args.results)
    comparison=json.loads((source/'comparison.json').read_text())
    evidence=json.loads((source/'finalization.evidence.json').read_text())
    state=validate_choices(json.loads(Path(args.decisions or source/'finalization.decisions.json').read_text()),evidence['spec'])
    if manifest['source_review_id'] != evidence['spec']['review_id']: raise ValueError('Source review identity differs')
    pending=sum(not c['reviewed'] for c in state['choices'].values())
    if pending and not args.preview:
        raise ValueError('Final primary assignments are not reviewed. Supply exported final decisions or use --preview explicitly.')
    models,proteins,coding=load_models(source,comparison,evidence)
    primary,cases=rows_for(evidence['spec'],models,state)
    extra=reference_panel(args.reference_panel,comparison['groups'])
    target=Path(args.out).resolve()
    if target==source or target in source.parents or source in target.parents:
        raise ValueError('Use a separate STOP 4 output directory')
    root=new_output(target); runner=TimedRunner(root)
    by_id={m['model_id']:m for m in models}; references=read_table(source/'reference_comparisons.tsv')
    synteny=read_table(source/'synteny.tsv'); rows=[]
    for row in primary+cases:
        result=assessment(by_id[row['model_id']],references,synteny)
        rows.append(dict(row,**result,protein_sha256=by_id[row['model_id']]['protein_sha256'],
                         cds_sha256=by_id[row['model_id']]['cds_sha256']))
    (root/'sequences').mkdir()
    for r in rows:
        mid=r['model_id']
        r['protein_file']='sequences/'+mid+'.protein.faa';r['cds_file']='sequences/'+mid+'.cds.fna'
        write_fasta(root/r['protein_file'],[(mid,proteins[mid])])
        write_fasta(root/r['cds_file'],[(mid,coding[mid])])
    for prefix,part in [('primary_targets',primary),('special_cases',cases)]:
        write_fasta(root/(prefix+'.proteins.faa'),[(r['model_id'],proteins[r['model_id']]) for r in part])
        write_fasta(root/(prefix+'.cds.fna'),[(r['model_id'],coding[r['model_id']]) for r in part])
    dump(root/'finalization.decisions.json',state)
    sequences=dict(proteins); sequence_groups={m['model_id']:m['group'] for m in models}; trees=[]
    for group in comparison['groups']:
        gid=group['id']; folder=root/'trees'/gid; folder.mkdir(parents=True)
        aligned=read_fasta(source/'queries'/gid/'proteins.aligned.faa')
        chosen={m['model_id'] for m in models if m['group']==gid}; qids={'reference_'+q for q in group['queries']}
        if not (chosen | qids) <= set(aligned): raise ValueError('Group MSA is missing a selected model or query')
        if any(aligned[mid].replace('-','') != proteins[mid] for mid in chosen):
            raise ValueError('Stored MSA differs from the exact selected protein')
        aligned={k:v for k,v in aligned.items() if k in chosen|qids}
        labels={mid:by_id[mid]['species']+' · '+by_id[mid]['original_model'] for mid in chosen}
        for qid in sorted(qids):
            seq=aligned[qid].replace('-','');sequences[qid]=seq
            if qid in sequence_groups and sequence_groups[qid]!=gid: raise ValueError('A reference occurs in multiple comparison groups; use distinct reference IDs')
            sequence_groups[qid]=gid
            info=next((r for r in references if r['reference_id']==qid and r['kind']=='protein'),{})
            labels[qid]=info.get('reference_species','Query')+' · '+info.get('reference_accession',qid)+' [query]'
        additions=[r for r in extra if r['group']==gid]
        for r in additions:
            sid=r['id']
            if sid in sequences: raise ValueError('Extra reference collides with saved sequence: '+sid)
            sequences[sid]=r['sequence'];sequence_groups[sid]=gid
            labels[sid]=r.get('species','Reference')+' · '+r.get('gene',sid)+' ['+r['role']+']'
        if additions:
            write_fasta(folder/'unaligned.faa',[(k,v.replace('-','')) for k,v in aligned.items()]+[(r['id'],r['sequence']) for r in additions])
            runner.run(['mafft','--auto','--thread',str(args.threads),folder/'unaligned.faa'],folder/'proteins.aligned.faa')
            aligned=read_fasta(folder/'proteins.aligned.faa')
            if set(aligned)!=set(labels) or any(aligned[k].replace('-','') != sequences[k] for k in aligned): raise ValueError('Extra-reference alignment changed input proteins')
        else: write_fasta(folder/'proteins.aligned.faa',aligned.items())
        outer=[r['id'] for r in additions if r['role']=='outgroup']
        if args.tree=='none': tree=dict(status='not_requested',reason='Tree calculation explicitly disabled.',leaves=[])
        else: tree=build_tree(folder,aligned,labels,runner,args.threads,args.bootstrap,args.occupancy,outer,
                             '; '.join(r['root_basis'] for r in additions if r['role']=='outgroup'),args.iqtree)
        tree.update(group=gid,gene=group['name'],folder=str(folder.relative_to(root)),
                    competitor_references=sum(r['role']=='paralog' for r in additions),
                    history_status='not_reconstructed',sequence_roles={k:('selected_model' if k in chosen else next((r['role'] for r in additions if r['id']==k),'supplied_query')) for k in aligned})
        dump(folder/'evidence.json',tree);trees.append(tree)
    catalogue=args.structures;retrieval=[]
    if getattr(args,'structure_provider','none')=='afdb':
        from .alphafold_db import collect
        # Chosen primary models get the pilot budget first, then separate
        # cases and supplied references; identical proteins use one request.
        order=list(dict.fromkeys([r['model_id'] for r in primary+cases]+list(sequences)))
        catalogue,retrieval=collect({sid:sequences[sid] for sid in order},root/'afdb',
            args.structure_cache,args.structure_limit,args.structure_offline,args.structure_refresh)
    elif getattr(args,'structure_provider',None)=='af3-prepare':
        from .af3_requests import prepare
        order=list(dict.fromkeys([r['model_id'] for r in primary+cases]+list(sequences)))
        retrieval=prepare({sid:sequences[sid] for sid in order},root/'af3',args.structure_limit)
    imported,structure_notes=import_structures(catalogue,sequences,root)
    retrieval_by_id={r['sequence_id']:r for r in retrieval}
    structural=compare_structures(imported,sequence_groups,root,runner,args.usalign)
    for r in rows:
        mid=r['model_id'];r['structure']=imported.get(mid)
        r['structure_retrieval']=retrieval_by_id.get(mid)
        possible=[]
        for p in structural:
            if p['status']!='computed' or p['scope']!='plddt_ge70' or p['same_protein_sequence'] or mid not in [p['sequence_a'],p['sequence_b']]:continue
            other=p['sequence_b'] if p['sequence_a']==mid else p['sequence_a']
            if other in by_id:continue  # Main row compares with a reference, not whichever candidate happens to fit best.
            side='a' if p['sequence_a']==mid else 'b'; refside='b' if side=='a' else 'a'
            possible.append(dict(reference_id=other,reference_label=next(t for t in trees if t['group']==r['group']).get('sequence_roles',{}).get(other,'reference'),
                tm_model=p['tm_score_'+side],tm_reference=p['tm_score_'+refside],coverage_model=p['coverage_full_'+side],coverage_reference=p['coverage_full_'+refside],rmsd_angstrom=p['rmsd_angstrom'],output_file=p['output_file']))
        r['structure_match']=max(possible,key=lambda v:min(v['tm_model'],v['tm_reference']),default=None)
    flat=[]
    for r in rows:
        v={k:v for k,v in r.items() if not isinstance(v,(dict,list))}
        v.update(facts=' '.join(r['facts']),limitations=' '.join(r['limitations']),
            structure_status='Exact structure available' if r['structure'] else (r['structure_retrieval']['status'] if r['structure_retrieval'] else 'Not assessed'),
            structure_missing_reason=r['structure_retrieval']['reason'] if r['structure_retrieval'] else '',
            structure_source=r['structure']['source'] if r['structure'] else '',
            structure_comparison_status='Computed on confident residues' if r['structure_match'] else 'No independent confident reference comparison')
        for k in ['tm_model','tm_reference','coverage_model','coverage_reference','rmsd_angstrom']:
            v['structure_'+k]=r['structure_match'][k] if r['structure_match'] else ''
        flat.append(v)
    for name,classification in [('final_table.tsv','primary_target'),('special_cases.tsv','special_case')]:
        part=[r for r in flat if r['classification']==classification]
        write_tsv(root/name,part,list(flat[0]) if flat else ['gene','species','model_id'])
    payload=dict(schema='stage4_dossier_v1',prototype=True,preview=bool(pending),
        review_id=evidence['spec']['review_id'],set_id=evidence['spec']['set_id'],
        rows=rows,trees=trees,structures=imported,structure_comparisons=structural,
        structure_notes=structure_notes,structure_retrieval=retrieval,counts=dict(primary=len(primary),special_cases=len(cases),
            species=len({r['genome_id'] for r in primary}),unreviewed_primary_buckets=pending),
        limitations=['No automatic orthology or functional-equivalence claim.',
            'Trees do not automatically validate target-versus-paralog identity or reconstruct duplication/loss history.',
            'No new folding jobs are submitted. Optional AFDB retrieval supplies existing full exact-sequence structures; absent exact structures require separate prediction.',
            'Confidence-filtered structural comparisons are not domain-specific comparisons. Low-confidence/disordered sequence is not diagnosed as a model error.'])
    dump(root/'evidence.json',payload)
    if imported:
        engine=args.usalign or shutil.which('USalign')
        if engine and any(c['status']=='computed' and c['scope']=='plddt_ge70' for c in structural):
            from .stage4_similarity import build_similarity
            payload['structure_geometry']=build_similarity(root,payload,engine)
        from .stage4_viewer import write_viewer,vendor_files
        if all(p.is_file() for p in vendor_files()):
            payload['structure_viewer']=write_viewer(root,payload)
            dump(root/'evidence.json',payload)
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<','\\u003c').replace('&','\\u0026')
    (root/'report.html').write_text(Path(__file__).with_name('stage4_report.txt').read_text().replace('__DATA__',text))
    (root/'START_HERE.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>STOP 4</title><p><a href="report.html">Open STOP 4 candidate dossier</a></p><p>Extract the whole folder. Opens offline in a modern browser.</p></html>')
    dump(root/'manifest.json',dict(stage='stage4_prototype',status='completed',created=now(),
        biological_approval_asserted=False,preview=bool(pending),source_manifest_sha256=sha256(source/'manifest.json'),
        source_seal_sha256=sha256(source/'seal.json'),final_decisions_sha256=sha256(root/'finalization.decisions.json'),
        input_catalogue_sha256=sha256(source/'sequence_catalogue.tsv'),
        extra_reference_panel_sha256=sha256(args.reference_panel) if args.reference_panel else None,
        structures_catalogue_sha256=sha256(args.structures) if args.structures else None,
        counts=payload['counts'],elapsed_seconds=time.perf_counter()-started,commands=runner.commands,
        source_code_sha256={n:sha256(Path(__file__).with_name(n)) if n.endswith('.txt') else sha256(Path(importlib.import_module('genecontext.'+n[:-3]).__file__)) for n in ['stage4.py','alphafold_db.py','af3_requests.py','stage4_phylogeny.py','stage4_structure.py','stage4_report.txt','stage4_viewer.py','stage4_viewer.txt','stage4_similarity.py','stage3_finalization.py','sequence.py','tools.py','cli.py']},
        parameters=dict(tree=args.tree,bootstrap=args.bootstrap,occupancy=args.occupancy,threads=args.threads,
            seed=42,structure_plddt_threshold=70,structure_provider=getattr(args,'structure_provider',None) or 'none',
            structure_limit=getattr(args,'structure_limit',3)),limitations=payload['limitations']))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    print(json.dumps(dict(report=str(root/'report.html'),preview=bool(pending),counts=payload['counts'],elapsed_seconds=time.perf_counter()-started)))
    return 0
