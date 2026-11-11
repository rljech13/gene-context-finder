"""One editable target representative per gene/genome, with traceable exceptions."""
from collections import defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from .gene_identity import identity as reference_identity, model_identities, export_identities
from .sequence import read_fasta, write_fasta, write_tsv
from .stage3_analysis import read_table


SCHEMA='stage3_final_choices_v1'


def normalized(value):
    return str(value or '').strip().casefold()


def validate_choices(state, spec):
    if not isinstance(state,dict) or state.get('schema')!=SCHEMA:
        raise ValueError('Expected stage-3 final choices JSON')
    if state.get('review_id')!=spec['review_id'] or state.get('set_id')!=spec['set_id']:
        raise ValueError('Final choices belong to another reviewed sequence set')
    buckets={b['key']:b for b in spec['buckets']}
    choices=state.get('choices')
    if not isinstance(choices,dict) or set(choices)!=set(buckets):
        raise ValueError('Final choices must cover the exact gene/genome buckets')
    for key,choice in choices.items():
        if not isinstance(choice,dict) or choice.get('primary') not in ['',*buckets[key]['model_ids']]:
            raise ValueError('A primary model must belong to its own gene/genome bucket')
        if not isinstance(choice.get('reason'),str) or not choice['reason'].strip():
            raise ValueError('Supply the basis for every primary assignment, including an unresolved assignment')
        if type(choice.get('reviewed')) is not bool:
            raise ValueError('The reviewed flag must be a boolean')
    notes=state.get('cases',{})
    allowed={m for b in buckets.values() for m in b['model_ids']}
    if not isinstance(notes,dict) or set(notes)-allowed:
        raise ValueError('Unknown model in special-case notes')
    for note in notes.values():
        if not isinstance(note,dict) or not isinstance(note.get('reason'),str) or not note['reason'].strip() or type(note.get('reviewed')) is not bool:
            raise ValueError('Special cases require a reason and boolean reviewed flag')
    return deepcopy(state)


def case_reason(model, primary, models):
    a=model['assignment']; identity=model['copy_identity']
    if a['known_other_copy']:
        return 'Specific source label is '+identity['display_name']+'; it is not automatically renamed to '+model['target_gene']+'.'
    if a['alternate_preferred']:
        ref=a['best_alternate']
        return 'Reference preference favors '+ref['symbol']+' ('+ref['accession']+'); kept separately as a possible related copy, not assigned cross-species orthology.'
    if not primary:
        return 'No target representative assigned in this organism; copy identity remains a separate case.'
    return 'Additional reviewed copy; '+models[primary]['original_model']+' is the current target representative in this organism.'


def rows_for(spec, models, state):
    by_id={m['model_id']:m for m in models};primary=[];cases=[]
    for bucket in spec['buckets']:
        choice=state['choices'][bucket['key']]
        for mid in bucket['model_ids']:
            m=by_id[mid];identity=m['copy_identity'];is_primary=choice['primary']==mid
            note=choice if is_primary else state.get('cases',{}).get(mid)
            reason=note['reason'] if note else case_reason(m,choice['primary'],by_id)
            ref=m['assignment']['best_target'] or {}
            row=dict(gene=bucket['gene'],group=bucket['group'],species=m['species'],genome_id=m['genome_id'],
                assembly_accession=m.get('assembly_accession',''),model_id=mid,original_model=m['original_model'],
                original_gene=m.get('original_gene',''),copy_label=m.get('copy_label',''),
                copy_identity=identity['display_name'],copy_identity_status=identity['status'],copy_identity_basis=identity['basis'],
                ncbi_gene_id=identity['ncbi_gene_id'],source=m['source'],contig=m['contig'],start=m['start'],end=m['end'],strand=m['strand'],
                protein_aa=m['protein_aa'],classification='primary_target' if is_primary else 'special_case',
                selection_reason=reason,decision_origin='reviewed' if note and note['reviewed'] else 'suggestion',
                best_target_reference=ref.get('accession',''),reference_gene=ref.get('symbol',''),
                protein_identity=ref.get('identity',''),model_coverage=ref.get('coverage_model',''),reference_coverage=ref.get('coverage_reference',''),
                evidence=' '.join(m['assignment']['facts']),cds_file=m.get('cds_file',''),protein_file=m.get('protein_file',''))
            (primary if is_primary else cases).append(row)
    return primary,cases


def prepare(root,payload,models,references):
    root=Path(root);data=json.loads((root/'source.evidence.json').read_text())
    identities=model_identities(data)
    export_identities(root,data['models'],identities,'model_id')
    group_names={g['id']:g['name'] for g in payload['groups']}
    for m in models:
        m['target_gene']=group_names[m['group']];m['copy_identity']=identities[m['model_id']]
    for group in payload['groups']:
        gid=group['id'];gene=group['name'];panel=references[gid]['protein']
        for ref in panel:
            member=(ref.get('aliases') or [{}])[0]
            ref['copy_identity']=reference_identity(dict(member,source='annotation',annotation_name=member.get('gene_symbol','')))
            ref['display_symbol']=ref['copy_identity']['label'] or ref['symbol']
        canonical={r['id'] for r in panel if normalized(r['symbol'])==normalized(gene)}
        fallback=not canonical
        if fallback:canonical={r['id'] for r in panel}
        pairs=read_table(root/'queries'/gid/'protein_pairs.tsv')
        for m in [m for m in models if m['group']==gid]:
            hits=[]
            for ref in panel:
                pair=next((p for p in pairs if {p['sequence_a'],p['sequence_b']}=={m['model_id'],ref['id']}),None)
                if not pair or pair['identity']=='':continue
                forward=pair['sequence_a']==m['model_id']
                cm=float(pair['coverage_a'] if forward else pair['coverage_b']);cr=float(pair['coverage_b'] if forward else pair['coverage_a'])
                hits.append(dict(reference_id=ref['id'],accession=ref['accession'],species=ref['species'],symbol=ref.get('display_symbol',ref['symbol']),
                    identity=float(pair['identity']),coverage_model=cm,coverage_reference=cr,score=float(pair['identity'])*min(cm,cr)))
            target=max((h for h in hits if h['reference_id'] in canonical),key=lambda h:h['score'],default=None)
            alternate=max((h for h in hits if h['reference_id'] not in canonical),key=lambda h:h['score'],default=None)
            identity=m['copy_identity'];known=identity['status'] in {'curated_locus_label','annotation_label','annotation_structure_match'}
            other=known and not fallback and normalized(identity['label'])!=normalized(gene)
            complete=m['quality']['coding_status']=='complete_predicted' and str(m.get('window_edge')).lower()!='true'
            facts=['Copy label: '+identity['display_name']+' ('+identity['basis'].replace('_',' ')+').',
                'Complete coding model.' if complete else 'Coding completeness or saved-window boundary needs review.']
            if target:facts.append('Best '+('supplied' if fallback else gene)+' query: '+target['accession']+'; identity '+format(target['identity'],'.1%')+', coverage model/reference '+format(target['coverage_model'],'.1%')+' / '+format(target['coverage_reference'],'.1%')+'.')
            else:facts.append('No comparable target-reference protein measurement.')
            agreement=m['quality']['source_agreement']
            facts.append({'exact_CDS_match':'Annotation and Tiberius have an identical CDS.',
                          'different_CDS_models':'Annotation and Tiberius propose different CDS models; the STOP-2 selection is retained.',
                          'other_source_unavailable':'Only one source is available for this model.'}.get(agreement,agreement))
            if identity['status'] in {'homology_only','unresolved'}:facts.append('Specific copy identity is not established by the source label.')
            if payload.get('reviewed_pair_count',0)<2:facts.append('Confirmed synteny was not used to establish this copy identity.')
            m['assignment']=dict(best_target=target,best_alternate=alternate,known_other_copy=other,complete=complete,
                alternate_preferred=bool(target and alternate and alternate['score']>target['score']),facts=facts,
                ranking=[int(complete),target['score'] if target else -1,int(known and not other),int(agreement=='exact_CDS_match')])
    groups=defaultdict(list)
    for m in models:groups[m['group'],m['genome_id']].append(m)
    buckets=[];choices={}
    for (group,genome),candidates in groups.items():
        key='bucket_'+hashlib.sha256((group+'\t'+genome).encode()).hexdigest()[:20]
        eligible=[m for m in candidates if not m['assignment']['known_other_copy'] and m['assignment']['best_target']]
        best=max(eligible,key=lambda m:m['assignment']['ranking'],default=None)
        if best:
            ref=best['assignment']['best_target']
            reason='Suggested representative: '+ref['accession']+' gives '+format(ref['identity'],'.1%')+' protein identity with model/reference coverage '+format(ref['coverage_model'],'.1%')+' / '+format(ref['coverage_reference'],'.1%')+'. '
            reason+=('Ranked first by coding completeness, then target-query identity × minimum coverage among the selected copies in this organism.' if len(candidates)>1 else 'Only STOP-2 selected copy in this organism.')
        else:reason='No primary target proposed: available copies have a different specific source label or lack target-reference comparison.'
        buckets.append(dict(key=key,group=group,gene=group_names[group],genome_id=genome,species=candidates[0]['species'],model_ids=[m['model_id'] for m in candidates]))
        choices[key]=dict(primary=best['model_id'] if best else '',reason=reason,reviewed=False)
    signature=[(m['group'],m['genome_id'],m['model_id'],m['cds_sha256']) for m in models]
    set_id=hashlib.sha256(json.dumps([payload['review_id'],sorted(signature)],sort_keys=True).encode()).hexdigest()
    spec=dict(review_id=payload['review_id'],set_id=set_id,buckets=buckets,
        method='One proposed representative per comparison group/genome. Specific known non-target labels are kept separate. Complete models rank first, then target-query identity × minimum model/reference coverage. Source label and exact CDS agreement break ties. This is a proposal, not resolved orthology.')
    state=dict(schema=SCHEMA,review_id=payload['review_id'],set_id=set_id,choices=choices,cases={})
    path=root/'finalization.decisions.json'
    if path.exists():state=validate_choices(json.loads(path.read_text()),spec)
    state=validate_choices(state,spec)
    primary,cases=rows_for(spec,models,state)
    fields=list((primary or cases)[0]) if primary or cases else ['gene','species','model_id','selection_reason']
    for name,rows in [('primary_targets.tsv',primary),('special_cases.tsv',cases),('final_model_choices.tsv',primary)]:write_tsv(root/name,rows,fields)
    raw={kind:read_fasta(root/name) if (root/name).stat().st_size else {} for kind,name in [('protein','selected.proteins.faa'),('cds','selected.cds.fna')]}
    for name,rows in [('primary_targets',primary),('special_cases',cases)]:
        for kind,suffix in [('protein','proteins.faa'),('cds','cds.fna')]:write_fasta(root/(name+'.'+suffix),[(r['model_id'],raw[kind][r['model_id']]) for r in rows])
    path.write_text(json.dumps(state,indent=2,ensure_ascii=False)+'\n')
    (root/'finalization.evidence.json').write_text(json.dumps(dict(spec=spec,models={m['model_id']:dict(copy_identity=m['copy_identity'],assignment=m['assignment']) for m in models}),indent=2,ensure_ascii=False)+'\n')
    return dict(spec=spec,state=state,rows=dict(primary=primary,cases=cases),raw_sequences=raw)
