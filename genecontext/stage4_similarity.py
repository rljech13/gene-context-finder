"""Per-residue deviation in the actual US-align superposition, not pLDDT."""
import json
import math
from pathlib import Path
import re
import subprocess
from .sequence import sha256
from .stage4_structure import ca_atoms,parse_usalign


def matrix_from_text(text):
    rows=[list(map(float,m.groups()[1:])) for line in text.splitlines()
          if (m:=re.fullmatch(r'\s*([012])\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*',line))]
    if len(rows)!=3 or not all(math.isfinite(v) for row in rows for v in row):raise ValueError('Invalid US-align transform')
    # Reflections and non-rigid transforms must never enter a similarity scene.
    u=[r[1:] for r in rows]
    for i in range(3):
        for j in range(3):
            if abs(sum(u[i][k]*u[j][k] for k in range(3))-(i==j))>1e-6:raise ValueError('Non-orthogonal US-align transform')
    det=sum(u[0][k]*(u[1][(k+1)%3]*u[2][(k+2)%3]-u[1][(k+2)%3]*u[2][(k+1)%3]) for k in range(3))
    if abs(det-1)>1e-6:raise ValueError('US-align reflection is not supported')
    return rows


def xyz(atom):return [float(atom['line'][a:b]) for a,b in [(30,38),(38,46),(46,54)]]
def move(point,matrix):return [row[0]+sum(row[k+1]*point[k] for k in range(3)) for row in matrix]
def key(atom):
    line=atom['line'];return line[21:22]+'|'+str(int(line[22:26]))+line[26:27].strip()


def residue_distances(text,matrix,aa,bb):
    lines=text.splitlines();i=next((i for i,l in enumerate(lines) if l.startswith('(":" denotes residue pairs')),None)
    if i is None or i+3>=len(lines):raise ValueError('Missing US-align residue correspondence')
    a,marks,b=lines[i+1:i+4]
    if len(a)!=len(b) or len(a)!=len(marks):raise ValueError('Invalid correspondence columns')
    if a.replace('-','')!=''.join(x['aa'] for x in aa) or b.replace('-','')!=''.join(x['aa'] for x in bb):raise ValueError('US-align correspondence sequence differs from input')
    ia=ib=0;pairs=[]
    for x,mark,y in zip(a,marks,b):
        if x!='-':ia+=1
        if y!='-':ib+=1
        if x!='-' and y!='-':
            if mark not in ':.':raise ValueError('Unmarked structural residue pair')
            distance=math.dist(move(xyz(aa[ia-1]),matrix),xyz(bb[ib-1]))
            pairs.append(dict(a=key(aa[ia-1]),b=key(bb[ib-1]),distance_angstrom=distance))
    stats=parse_usalign(text)
    if len(pairs)!=stats['aligned_residues']:raise ValueError('Structural pair count differs')
    # Independent check against US-align's -do distances (rounded to 0.001 A).
    output=[];collect=False
    for line in lines:
        if line.startswith('#Aligned atom 1'):collect=True;continue
        if collect and line.startswith('###'):break
        if collect and line.strip():
            fields=line.split('\t')
            if len(fields)!=3:raise ValueError('Malformed US-align distance table')
            output.append(float(fields[-1]))
    if len(output)!=len(pairs) or any(abs(p['distance_angstrom']-d)>.0011 for p,d in zip(pairs,output)):
        raise ValueError('Displayed distances do not reproduce US-align distance table')
    return pairs,stats


def build_similarity(root,data,engine):
    """Reproduce saved unique confident-region pairs with explicit matrices."""
    root=Path(root);structures=data['structures'];selected={r['model_id'] for r in data['rows']}
    representatives={}
    for sid,s in sorted(structures.items(),key=lambda item:item[0] not in selected):representatives.setdefault(s['protein_sha256'],sid)
    aliases={sid:representatives[s['protein_sha256']] for sid,s in structures.items()}
    comparisons=[];seen=set();folder=root/'structure_similarity';folder.mkdir(exist_ok=True)
    for old in data.get('structure_comparisons',[]):
        if old['status']!='computed' or old['scope']!='plddt_ge70':continue
        a,b=aliases[old['sequence_a']],aliases[old['sequence_b']]
        pair=tuple(sorted([a,b]))
        if a==b or pair in seen:continue
        seen.add(pair)
        pa,pb=[root/structures[s]['confident_file'] for s in [a,b]]
        if not pa.is_file() or not pb.is_file():raise ValueError('Missing confident structure for saved pair')
        name=a+'.'+b;matrix=folder/(name+'.matrix.txt');raw=folder/(name+'.usalign.txt')
        command=[str(engine),str(pa),str(pb),'-m',str(matrix),'-do']
        run=subprocess.run(command,capture_output=True,text=True,check=True,timeout=120)
        raw.write_text(run.stdout);(folder/(name+'.stderr.txt')).write_text(run.stderr)
        transform=matrix_from_text(matrix.read_text())
        pairs,stats=residue_distances(run.stdout,transform,ca_atoms(pa),ca_atoms(pb))
        for field in ['aligned_residues','rmsd_angstrom','tm_score_a','tm_score_b']:
            if abs(stats[field]-old[field])>1e-6:raise ValueError('Saved structural measurement changed: '+field)
        entry=dict(target=a,reference=b,scope='plddt_ge70',matrix=transform,pairs=pairs,**stats,
            coverage_full_a=len(pairs)/structures[a]['full_length'],coverage_full_b=len(pairs)/structures[b]['full_length'],
            display_rmsd_angstrom=math.sqrt(sum(p['distance_angstrom']**2 for p in pairs)/len(pairs)),
            raw_file=str(raw.relative_to(root)),matrix_file=str(matrix.relative_to(root)),
            command=command,input_sha256=[sha256(pa),sha256(pb)],engine_sha256=sha256(engine))
        comparisons.append(entry)
    data['structure_similarity']=comparisons
    (folder/'comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    return dict(pairs=len(comparisons),scope='plddt_ge70',engine_sha256=sha256(engine),distance_validation='All C-alpha distances reproduce US-align -do within 0.0011 A')
