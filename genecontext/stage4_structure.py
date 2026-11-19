"""Exact-sequence structure imports and measured, optional structural comparisons."""
from itertools import combinations
import json
import math
from pathlib import Path
import re
import shutil

from .sequence import sha256, write_tsv

RESIDUES = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))


def ca_atoms(path, chain='A'):
    atoms, seen = [], set()
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA' or line[21:22] != chain: continue
        if line[16:17] not in (' ', 'A'): continue
        key = line[22:27]
        if key in seen: raise ValueError('Duplicate CA residue in structure')
        seen.add(key)
        residue = line[17:20].strip()
        if residue not in RESIDUES: raise ValueError('Unsupported residue in prototype PDB: '+residue)
        xyz = [float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]
        confidence = float(line[60:66])
        if not all(math.isfinite(x) for x in [*xyz,confidence]): raise ValueError('Non-finite PDB coordinate/confidence')
        atoms.append(dict(aa=RESIDUES[residue], line=line, confidence=confidence))
    if not atoms: raise ValueError('No CA atoms in requested structure chain')
    return atoms


def import_structures(catalogue, sequences, root):
    if not catalogue: return {}, []
    file = Path(catalogue).resolve(); data = json.loads(file.read_text())
    if data.get('schema') != 'stage4_structures_v1': raise ValueError('Expected stage4_structures_v1 catalogue')
    imported, notes = {}, []
    for entry in data['entries']:
        sid = entry['sequence_id']
        if sid not in sequences or sid in imported: raise ValueError('Unknown or duplicate structure sequence ID: '+sid)
        if not entry.get('source'): raise ValueError('Structure provenance source is required')
        path = (file.parent/entry['file']).resolve()
        if sha256(path) != entry['file_sha256']: raise ValueError('Structure file hash mismatch: '+sid)
        import hashlib
        if hashlib.sha256(sequences[sid].encode()).hexdigest() != entry['protein_sha256']:
            raise ValueError('Structure catalogue protein hash mismatch: '+sid)
        atoms = ca_atoms(path, entry.get('chain','A'))
        if ''.join(a['aa'] for a in atoms) != sequences[sid]:
            raise ValueError('Structure CA sequence is not the exact full selected protein: '+sid)
        confidence = entry.get('confidence') == 'plddt_bfactor'
        if confidence and any(not 0 <= a['confidence'] <= 100 for a in atoms): raise ValueError('Invalid pLDDT')
        dest = root/'structures'/(sid+'.pdb'); dest.parent.mkdir(exist_ok=True)
        dest.write_text('\n'.join(a['line'] for a in atoms)+'\nEND\n')
        full = dest.with_name(sid+'.full.pdb')
        full_atoms = [line for line in path.read_text().splitlines()
                      if line.startswith('ATOM  ') and line[21:22]==entry.get('chain','A') and line[16:17] in (' ','A')]
        full.write_text('\n'.join(full_atoms)+'\nEND\n')
        confident = [a for a in atoms if a['confidence'] >= 70] if confidence else []
        masked = dest.with_name(sid+'.confident.pdb')
        if confident: masked.write_text('\n'.join(a['line'] for a in confident)+'\nEND\n')
        imported[sid] = dict(sequence_id=sid, source=entry['source'], accession=entry.get('accession',''),
            source_url=entry.get('source_url',''), protein_sha256=entry['protein_sha256'],
            input_file_sha256=entry['file_sha256'], file=str(dest.relative_to(root)),
            file_sha256=sha256(dest), full_length=len(atoms),
            full_file=str(full.relative_to(root)),full_file_sha256=sha256(full),
            coordinate_scope='all_atom_chain' if len(full_atoms)>len(atoms) else 'CA_trace',
            confidence_kind='pLDDT' if confidence else 'not_available',
            confident_residues=len(confident) if confidence else None,
            fraction_plddt_ge70=len(confident)/len(atoms) if confidence else None,
            mean_plddt=sum(a['confidence'] for a in atoms)/len(atoms) if confidence else None,
            confident_file=str(masked.relative_to(root)) if confident else '',
            pae_status='not_imported', exact_sequence_match=True)
    notes.append('PDB CA-chain imports require a full, exact sequence match. pLDDT is used only when explicitly declared; experimental B factors are not pLDDT. PAE and domain boundaries are not assessed in this prototype.')
    return imported, notes


def parse_usalign(text):
    aligned = re.search(r'Aligned length\s*=\s*(\d+),\s*RMSD\s*=\s*([\d.eE+-]+)', text)
    tm = re.findall(r'TM-score\s*=\s*([\d.eE+-]+)[^\n]*normalized by length of (?:Structure|Chain)_([12])', text)
    if not aligned or {n for _,n in tm} != {'1','2'}: raise ValueError('Incomplete US-align measurements')
    values = dict(aligned_residues=int(aligned[1]), rmsd_angstrom=float(aligned[2]),
                  tm_score_a=float(next(v for v,n in tm if n=='1')), tm_score_b=float(next(v for v,n in tm if n=='2')))
    if not all(math.isfinite(v) for v in values.values()) or values['rmsd_angstrom'] < 0 or not all(0<=values[k]<=1 for k in ['tm_score_a','tm_score_b']):
        raise ValueError('Invalid US-align measurements')
    return values


def compare_structures(imported, sequence_groups, root, runner, engine=None):
    engine = engine or shutil.which('USalign')
    rows = []
    for a,b in combinations(imported,2):
        if sequence_groups[a] != sequence_groups[b]: continue
        for scope in ['full_chain','plddt_ge70']:
            aa,bb = imported[a],imported[b]
            pa,pb = (aa['file'],bb['file']) if scope=='full_chain' else (aa['confident_file'],bb['confident_file'])
            na,nb = (aa['full_length'],bb['full_length']) if scope=='full_chain' else (aa['confident_residues'],bb['confident_residues'])
            row = dict(sequence_a=a,sequence_b=b,group=sequence_groups[a],scope=scope,
                status='not_assessed',reason='',tm_score_a='',tm_score_b='',rmsd_angstrom='',aligned_residues='',
                coverage_full_a='',coverage_full_b='',retained_a=na,retained_b=nb,
                same_protein_sequence=aa['protein_sha256']==bb['protein_sha256'],output_file='')
            if not pa or not pb or min(na or 0,nb or 0) < 20:
                row['reason']='Fewer than 20 comparable CA residues or confidence unavailable.'
            elif not engine: row['reason']='US-align not installed.'
            else:
                output = root/'structures'/(a+'.'+b+'.'+scope+'.txt')
                runner.run([engine,root/pa,root/pb],output)
                value = parse_usalign(output.read_text())
                if not 0 < value['aligned_residues'] <= min(na,nb): raise ValueError('US-align aligned length exceeds imported coordinates')
                row.update(value,status='computed',coverage_full_a=value['aligned_residues']/aa['full_length'],
                           coverage_full_b=value['aligned_residues']/bb['full_length'],output_file=str(output.relative_to(root)))
            rows.append(row)
    write_tsv(root/'structure_comparisons.tsv',rows,list(rows[0]) if rows else ['sequence_a','sequence_b','status'])
    return rows
