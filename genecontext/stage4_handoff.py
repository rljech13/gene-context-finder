"""Consume the exact alignment-first STOP 3 export without rebuilding its MSA or tree."""
import hashlib
import json
from pathlib import Path
import re
import shutil
from .stage4_phylogeny import parse_newick, tip_names, tree_svg
from .tree_views import unrooted_svg

from .sequence import sha256, translate, write_fasta, write_tsv
from .stage1 import dump, new_output, now
from .stage3 import TimedRunner
from .stage4_structure import import_structures, compare_structures


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def safe_id(value):
    return isinstance(value, str) and re.fullmatch(r'[\w.-]{1,150}', value, re.ASCII) and value not in {'.', '..'}


def read_handoff(path):
    source = Path(path)
    if source.stat().st_size > 64 * 1024 * 1024:
        raise ValueError('STOP 4 input exceeds 64 MiB')
    envelope = json.loads(source.read_bytes())
    payload = envelope['payload']
    if sha(payload) != envelope.get('sha256'):
        raise ValueError('STOP 4 input checksum differs')
    data = json.loads(payload)
    stage2, stage3 = data['stage2'], data['stage3']
    if (data['schema'] != 'gcf_stop4_input_v1' or stage2['schema'] != 'stage2_draft_v1'
            or stage2['review_id'] != data['source_review_id'] or stage3['reviewed'] is not True):
        raise ValueError('Expected a reviewed alignment-first STOP 3 export')
    members = {s['id']: s for s in data['members']}
    if len(members) != len(data['members']) or any(not safe_id(i) for i in members):
        raise ValueError('Invalid or duplicate sequence IDs')
    for mid, member in members.items():
        meta = member['metadata']
        protein, cds = member['protein'], member.get('cds') or ''
        if not isinstance(member['group'],str) or not member['group'].strip() or not re.fullmatch('[A-Z]+', protein):
            raise ValueError('Invalid protein or group: ' + mid)
        if not member['reference'] and (
            meta['model_id'] != mid or meta['protein_sha256'] != sha(protein)
            or meta.get('cds_sha256', '') != sha(cds) or not cds
            or not re.fullmatch('[ACGTRYSWKMBDHVN]+', cds) or len(cds) % 3
            or translate(cds).removesuffix('*') != protein
            or meta['start'] < 0 or meta['end'] <= meta['start'] or meta['strand'] not in ('+', '-')
        ):
            raise ValueError('Model identity, translation or genomic coordinates differ: ' + mid)
        if member['reference'] and cds:
            raise ValueError('The reference panel is protein-only')
    chosen = {mid for mid, choice in stage2['models'].items()
              if choice['decision'] == 'inspect' and choice['role'] in ('target', 'alignment_candidate')}
    choices, analyses = stage3['decisions'], data['analyses']
    if not chosen or chosen != set(choices) or not chosen <= members.keys():
        raise ValueError('STOP 3 must account for each confirmed candidate')
    regions = {r['region_id']: r for r in stage3['synteny']['regions']}
    if len(regions) != len(stage3['synteny']['regions']):
        raise ValueError('Duplicate synteny regions')
    if set(regions) != {members[mid]['metadata']['region_id'] for mid in chosen}:
        raise ValueError('Synteny regions do not match the reviewed candidates')
    reviewed_pairs = stage3['synteny']['reviewed_pairs']
    for evidence in reviewed_pairs.values():
        if evidence.get('kind') != 'neighbor_proposal':
            raise ValueError('Unexpected type of reviewed neighbor evidence')
    used_analyses = set()
    physical, primary = set(), set()
    for mid, choice in choices.items():
        member = members[mid]
        group, meta = member['group'], member['metadata']
        if group not in analyses or member['reference']:
            raise ValueError('Final decisions require a prepared candidate and its group analysis')
        used_analyses.add(group)
        region = regions[meta['region_id']]
        if (region['genome_id'] != meta['genome_id'] or region['contig'] != meta['contig']
                or stage2['regions'][meta['region_id']]['decision'] != 'inspect'
                or stage2['regions'][meta['region_id']]['saved'] is not True):
            raise ValueError('Candidate context differs: ' + mid)
        role = choice['role']
        if role not in {'primary_target', 'special_case', 'unresolved', 'excluded'}:
            raise ValueError('Final model choice requires an outcome: ' + mid)
        if not isinstance(choice.get('reason', ''), str):
            raise ValueError('Invalid final model comment: ' + mid)
        if role in {'unresolved', 'excluded'}:
            continue
        if meta.get('coding_status') == 'disrupted' or mid not in analyses[group]['alignment']['ids']:
            raise ValueError('Retained model is disrupted or absent from the saved alignment: ' + mid)
        locus = (meta['genome_id'], meta.get('model_constraints', {}).get('physical_id') or meta.get('gene_group') or meta['region_id'])
        if locus in physical:
            raise ValueError('Multiple final models for one physical gene')
        physical.add(locus)
        if role == 'primary_target':
            native = stage2['models'][mid]
            if (meta.get('model_constraints', {}).get('target_warnings')
                    and not (native.get('target_override') is True and native.get('notes','').strip())):
                raise ValueError('Primary model has an unresolved target warning: ' + mid)
            bucket = (group, meta['genome_id'])
            if bucket in primary or mid not in analyses[group]['tree']['ids']:
                raise ValueError('A primary model must be unique in the genome/group and included in the tree input')
            primary.add(bucket)
    if used_analyses != set(analyses):
        raise ValueError('Unexpected comparison group')
    rows_used = set()
    for group, part in analyses.items():
        alignment, tree = part['alignment'], part['tree']
        comments = part.get('comments', {})
        if (not isinstance(comments, dict) or not comments.keys() <= {'alignment', 'tree'}
                or any(not isinstance(value, str) for value in comments.values())):
            raise ValueError('Invalid alignment or tree review comment: ' + group)
        if tree['alignment_id'] != alignment['id'] or not set(tree['ids']) <= set(alignment['ids']):
            raise ValueError('Tree provenance differs from its saved source alignment')
        for phase, checkpoint in [('alignment', alignment), ('tree', tree)]:
            ids = checkpoint['ids']
            if (phase not in {'alignment', 'tree'} or checkpoint['phase'] != phase
                    or checkpoint['group'] != group or checkpoint['review_id'] != data['alignment_review_id']
                    or not re.fullmatch(r'[0-9a-f]{20,64}', checkpoint['id'])
                    or len(ids) < 2 or len(ids) != len(set(ids)) or set(ids) != set(checkpoint['protein'])
                    or len({len(row) for row in checkpoint['protein'].values()}) != 1
                    or set(checkpoint['cds']) != {mid for mid in ids if members[mid].get('cds')}):
                raise ValueError('Alignment or tree review identity differs')
            rows_used.update(ids)
            for mid in ids:
                row = checkpoint['protein'][mid]
                if members[mid]['group'] != group or row.replace('-', '') != members[mid]['protein']:
                    raise ValueError('Protein sequence changed in saved alignment: ' + mid)
                if phase == 'tree' and (row != alignment['protein'][mid]
                        or checkpoint['cds'].get(mid) != alignment['cds'].get(mid)):
                    raise ValueError('Tree input columns differ from the saved alignment')
                if mid in checkpoint['cds']:
                    coding = checkpoint['cds'][mid]
                    native = members[mid]['cds']
                    if len(native) == (len(members[mid]['protein'])+1)*3:
                        native = native[:-3]
                    if (len(coding) != len(row)*3 or coding.replace('-', '') != native
                            or any((aa == '-') != (coding[i*3:i*3+3] == '---') for i, aa in enumerate(row))):
                        raise ValueError('CDS columns changed in the saved alignment: ' + mid)
        phy = tree['phylogeny']
        if phy['status'] not in {'computed', 'not_computed'} or (phy['status'] == 'computed' and not phy.get('unrooted_newick')):
            raise ValueError('Invalid saved tree status')
        if phy['status'] == 'not_computed' and not phy.get('reason'):
            raise ValueError('Missing reason for unavailable tree')
        leaves = phy['leaves']
        leaf_ids = [mid for leaf in leaves for mid in leaf['sequence_ids']]
        if len(leaf_ids) != len(set(leaf_ids)) or set(leaf_ids) != set(tree['ids']):
            raise ValueError('Tree aliases differ from the tree input')
        if any(sha(members[mid]['protein']) != leaf['protein_sha256'] for leaf in leaves for mid in leaf['sequence_ids']):
            raise ValueError('Tree alias protein digest differs')
        if phy.get('tree_input_fasta') is not None:
            mask = phy['mask']
            columns = mask['retained_columns_1based']
            width = len(next(iter(tree['protein'].values())))
            if (mask['original_columns'] != width or mask['retained_columns'] != len(columns)
                    or columns != sorted(set(columns)) or any(not 1 <= c <= width for c in columns)):
                raise ValueError('Tree column mask differs from the source alignment')
            alphabet = set('ACDEFGHIKLMNPQRSTVWY')
            expected_input = {}
            for leaf in leaves:
                row = tree['protein'][leaf['sequence_ids'][0]]
                projected = ''.join(row[c-1] if row[c-1] in alphabet else '-' for c in columns)
                if leaf['included'] != bool(set(projected)&alphabet):
                    raise ValueError('Tree leaf inclusion differs from retained columns')
                if leaf['included']:
                    expected_input[leaf['tree_id']] = projected
            fasta_rows = {}
            for entry in phy['tree_input_fasta'].split('>')[1:]:
                name, *seq = entry.splitlines()
                if name in fasta_rows:
                    raise ValueError('Duplicate tree-input FASTA row')
                fasta_rows[name] = ''.join(seq)
            if fasta_rows != expected_input:
                raise ValueError('Tree-input FASTA differs from the saved column mask')
        expected = {leaf['tree_id'] for leaf in leaves if leaf['included']}
        if len(expected) != sum(leaf['included'] for leaf in leaves) or any(not safe_id(t) for t in expected):
            raise ValueError('Invalid or duplicate tree leaf ID')
        for mode in ('unrooted', 'rooted'):
            text = phy.get(mode+'_newick')
            if text:
                tips = tip_names(parse_newick(text))
                if set(tips) != expected or len(tips) != len(expected):
                    raise ValueError('Newick tip set differs from the saved tree input')
    if rows_used | chosen != set(members):
        raise ValueError('Foreign sequence in handoff')
    return data, members


def run(args):
    data, members = read_handoff(args.handoff)
    if getattr(args,'structure_provider',None) not in (None,'none'):
        raise ValueError('Use structure-fetch separately and pass its exact --structures catalogue. Folding preparation requires separate post-synteny reference approval.')
    if args.decisions or args.reference_panel or args.preview or args.iqtree or args.tree != 'auto':
        raise ValueError('The STOP 3 export already fixes the decisions, alignment, tree and references')
    root = new_output(Path(args.out).resolve())
    (root/'sequences').mkdir()
    dump(root/'finalization.decisions.json', data['stage3']['decisions'])
    dump(root/'synteny.review.json', data['stage3']['synteny'])
    dump(root/'model.review.json', data['stage2'])
    trees = []
    for group, part in data['analyses'].items():
        folder = root/'trees'/('group_'+sha(group)[:16])
        folder.mkdir(parents=True)
        alignment, saved = part['alignment'], part['tree']
        write_fasta(folder/'proteins.aligned.faa', alignment['protein'].items())
        write_fasta(folder/'cds.codon.aligned.fna', alignment['cds'].items())
        phy = {k:v for k,v in saved['phylogeny'].items() if not k.endswith('_svg')}
        tree = dict(phy, group=group, gene=group, folder=str(folder.relative_to(root)),
                    alignment_result_id=alignment['id'], tree_result_id=saved['id'])
        if part.get('comments'):
            tree['review_comments'] = part['comments']
        if tree.get('tree_input_fasta'):
            (folder/'tree.aligned.faa').write_text(tree.pop('tree_input_fasta'))
            tree['used_alignment'] = 'tree.aligned.faa'
        for mode in ('unrooted', 'rooted'):
            if tree.get(mode+'_newick'):
                (folder/(mode+'.nwk')).write_text(tree[mode+'_newick'])
        if tree['status'] == 'computed':
            preferred = 'rooted' if tree.get('rooted_newick') else 'unrooted'
            tree['rooting'] = 'outgroup' if preferred == 'rooted' else 'unrooted'
            (folder/'ml.treefile').write_text(tree[preferred+'_newick'])
            # Render a new safe visualization of the preserved Newick, instead of copying active markup.
            (folder/'unrooted.svg').write_text(unrooted_svg(tree['unrooted_newick'], tree['leaves'], scale_lengths=True))
            if preferred == 'rooted':
                (folder/'rooted.svg').write_text(tree_svg(tree['rooted_newick'], tree['leaves']))
            shutil.copyfile(folder/(preferred+'.svg'), folder/'tree.svg')
        dump(folder/'evidence.json', tree)
        trees.append(tree)
    selected = {mid for mid, choice in data['stage3']['decisions'].items() if choice['role'] in {'primary_target','special_case'}}
    structural_pool = {mid:s for mid,s in members.items() if mid in selected or s['reference']}
    structures, notes = import_structures(args.structures, {mid:s['protein'] for mid,s in structural_pool.items()}, root)
    runner = TimedRunner(root)
    structural = compare_structures(structures, {mid:s['group'] for mid,s in structural_pool.items()}, root, runner,args.usalign)
    rows = []
    for mid, outcome in data['stage3']['decisions'].items():
        if mid not in selected:
            continue
        s = members[mid]
        m = s['metadata']
        protein_file, cds_file = 'sequences/'+mid+'.protein.faa', 'sequences/'+mid+'.cds.fna'
        write_fasta(root/protein_file, [(mid,s['protein'])])
        write_fasta(root/cds_file, [(mid,s['cds'])])
        hits = sorted(m.get('target_support') or [], key=lambda h:h['bitscore'], reverse=True)
        hit = hits[0] if hits else {}
        model_status = 'Coding consistent' if m['coding_status']=='complete_predicted' and not m.get('window_edge') else 'Check model'
        measurements = []
        for p in structural:
            if (p['status']!='computed' or p['scope']!='plddt_ge70'
                    or p['same_protein_sequence'] or mid not in [p['sequence_a'],p['sequence_b']]):
                continue
            side, ref = ('a','b') if p['sequence_a']==mid else ('b','a')
            other=p['sequence_'+ref]
            if other in members and members[other]['reference']:
                measurements.append(dict(reference_id=other, tm_model=p['tm_score_'+side],
                    tm_reference=p['tm_score_'+ref],coverage_model=p['coverage_full_'+side],
                    coverage_reference=p['coverage_full_'+ref],rmsd_angstrom=p['rmsd_angstrom'],output_file=p['output_file']))
        rows.append(dict(model_id=mid,group=s['group'],gene=s['group'],species=s['label'].split(' · ')[0],
            genome_id=m['genome_id'],original_model=m['original_model'],source=m['source'],contig=m['contig'],
            start=m['start'],end=m['end'],strand=m['strand'],decision_origin='reviewed',classification=outcome['role'],
            selection_reason=outcome.get('reason',''),synteny_status='Reviewed',copy_identity=m.get('annotation_name',''),
            identity_assessment='Reviewed candidate', model_assessment=model_status,model_notes=m['coding_status'],
            reference_accession=hit.get('subject'),protein_identity=hit.get('identity'),coverage_model=hit.get('coverage_query'),
            coverage_reference=hit.get('coverage_subject'),protein_aa=len(s['protein']),protein_file=protein_file,cds_file=cds_file,
            protein_sha256=m['protein_sha256'],cds_sha256=m['cds_sha256'],structure=structures.get(mid),
            structure_match=max(measurements,key=lambda p:min(p['tm_model'],p['tm_reference']),default=None),
            facts=['Exact selected protein and CDS verified.', 'Genome context reviewed in STOP 3.',
                   'Native coding status: '+m['coding_status']+'.'],
            
            limitations=['Tree, sequence and structure evidence are recorded separately. Orthology and expression require biological interpretation.']))
    counts={'unresolved':sum(v['role']=='unresolved' for v in data['stage3']['decisions'].values()),
            'excluded':sum(v['role']=='excluded' for v in data['stage3']['decisions'].values())}
    for name,role in [('primary_targets','primary_target'),('special_cases','special_case')]:
        part=[r for r in rows if r['classification']==role]
        counts[role]=len(part)
        write_fasta(root/(name+'.proteins.faa'), [(r['model_id'],members[r['model_id']]['protein']) for r in part])
        write_fasta(root/(name+'.cds.fna'), [(r['model_id'],members[r['model_id']]['cds']) for r in part])
        flat=[]
        for r in part:
            item={k:v for k,v in r.items() if not isinstance(v,(dict,list))}
            item.update(structure_source=(r['structure'] or {}).get('source',''),
                structure_reference=(r['structure_match'] or {}).get('reference_id',''))
            for k in ('tm_model','tm_reference','coverage_model','coverage_reference','rmsd_angstrom'):
                item['structure_'+k]=(r['structure_match'] or {}).get(k,'')
            flat.append(item)
        write_tsv(root/('final_table.tsv' if role=='primary_target' else 'special_cases.tsv'),flat,list(flat[0]) if flat else ['gene','species','model_id'])
    payload=dict(schema='stage4_dossier_v1',prototype=True,preview=False,review_id=data['source_review_id'],
        set_id=sha256(args.handoff),rows=rows,trees=trees,structures=structures,structure_comparisons=structural,
        structure_notes=notes,counts=counts)
    if structures:
        engine=args.usalign or shutil.which('USalign')
        if engine and any(c['status']=='computed' and c['scope']=='plddt_ge70' for c in structural):
            from .stage4_similarity import build_similarity
            payload['structure_geometry']=build_similarity(root,payload,engine)
        from .stage4_viewer import write_viewer, vendor_files
        if all(p.is_file() for p in vendor_files()):
            payload['structure_viewer']=write_viewer(root,payload)
    dump(root/'evidence.json',payload)
    text=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('<',chr(92)+'u003c').replace('&',chr(92)+'u0026')
    (root/'report.html').write_text(Path(__file__).with_name('stage4_report.txt').read_text().replace('__DATA__',text).replace('Evidence JSON</a> ·', 'Evidence JSON</a> · <a href="synteny.review.json">Reviewed synteny</a> ·'))
    (root/'START_HERE.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>STOP 4</title><p><a href="report.html">Open STOP 4 candidate dossier</a></p></html>')
    dump(root/'manifest.json',dict(stage='stage4_prototype',status='completed',created=now(),
        source_review_id=data['source_review_id'],alignment_review_id=data['alignment_review_id'],
        source_code_sha256={n:sha256(Path(__file__).with_name(n)) for n in ('stage4_handoff.py','stage4_report.txt','stage4_structure.py','stage4_similarity.py','stage4_viewer.py','tree_views.py')},
        biological_approval_asserted=False,source_handoff_sha256=sha256(args.handoff),counts=counts,commands=runner.commands))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file()})
    print(json.dumps(dict(report=str(root/'report.html'),counts=counts)))
    return 0
