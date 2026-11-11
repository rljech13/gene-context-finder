#!/usr/bin/env python3
"""Independent STOP-2 output audit. Run on host2 through submit.py.

Uses BioPython to reread genomic FASTA and translate native-coordinate CDS.
It does not select models, establish orthology, or modify scientific output.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def fasta(path):
    with open(path) as handle:
        records=list(SeqIO.parse(handle,'fasta'))
    assert len({r.id for r in records})==len(records),'Duplicate FASTA ID'
    return {r.id:str(r.seq).upper() for r in records}


def translation(dna):
    # The pipeline deliberately marks every ambiguous codon X. BioPython can
    # sometimes resolve an IUPAC codon more narrowly, so apply that policy here.
    return ''.join(str(Seq(dna[i:i+3]).translate()) if set(dna[i:i+3])<=set('ACGT') else 'X'
                   for i in range(0,len(dna)-2,3))


def native_cds(model,sequence):
    blocks=sorted(model['cds_blocks'],reverse=model['strand']=='-')
    def cut(a,b):
        assert 0<=a<b<=len(sequence),'Coding interval outside chromosome'
        value=sequence[a:b]
        return str(Seq(value).reverse_complement()) if model['strand']=='-' else value
    parts=[cut(a,b) for a,b,_ in blocks];parts[0]=parts[0][blocks[0][2] or 0:]
    joined=''.join(parts);extra=[]
    for start,end in model['stop_blocks']:
        runs=[(start,end)]
        for ca,cb,_ in blocks:
            pieces=[]
            for a,b in runs:
                if a<ca and min(b,ca)>a:pieces.append((a,min(b,ca)))
                if b>cb and max(a,cb)<b:pieces.append((max(a,cb),b))
            runs=pieces
        extra.extend(runs)
    if extra:
        tail=''.join(cut(a,b) for a,b in sorted(extra,reverse=model['strand']=='-'))
        terminal=min(a for a,b,_ in blocks) if model['strand']=='-' else max(b for a,b,_ in blocks)
        adjacent=max(b for a,b in extra)==terminal if model['strand']=='-' else min(a for a,b in extra)==terminal
        if len(tail)<=3 and adjacent and translation(joined+tail).endswith('*'):joined+=tail
    return joined


def raw_hits(path,queries,subjects):
    groups=defaultdict(list)
    for line in path.read_text().splitlines():
        f=line.split();assert len(f)==12
        q,s=f[:2];assert q in queries and s in subjects
        identity,length=float(f[2])/100,int(f[3]);qs,qe,ss,se=map(int,f[4:8]);e,score=map(float,f[8:10])
        assert all(math.isfinite(x) for x in [identity,e,score])
        assert 0<=identity<=1 and e>=0 and score>=0 and length>0
        assert (int(f[10]),int(f[11]))==(queries[q],subjects[s])
        assert 1<=qs<=qe<=queries[q] and 1<=ss<=se<=subjects[s]
        groups[q,s].append((qs-1,qe,ss-1,se,identity,length,e,score))
    result={}
    for key,hsps in groups.items():
        chosen=[]
        for h in sorted(hsps,key=lambda x:-x[7]):
            if any(max(h[0],a[0])<min(h[1],a[1]) or max(h[2],a[2])<min(h[3],a[3]) for a in chosen):continue
            chosen.append(h)
        q,s=key;n=sum(h[5] for h in chosen)
        result[key]=dict(bitscore=sum(h[7] for h in chosen),evalue=min(h[6] for h in chosen),
            identity=sum(h[4]*h[5] for h in chosen)/n,
            coverage_query=sum(h[1]-h[0] for h in chosen)/queries[q],
            coverage_subject=sum(h[3]-h[2] for h in chosen)/subjects[s])
    return result


def same_statistics(observed,expected):
    for key in ['bitscore','evalue','identity','coverage_query','coverage_subject']:
        assert math.isclose(observed[key],expected[key],rel_tol=1e-9,abs_tol=1e-12),(key,observed,expected)


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',required=True);p.add_argument('--source',required=True)
    p.add_argument('--draft',required=True);p.add_argument('--genomes',required=True);p.add_argument('--out',required=True)
    p.add_argument('--pilot');p.add_argument('--reuse-native-audit',help='Passed independent audit of the unchanged parent checkpoint')
    a=p.parse_args();root=Path(a.results);data=json.loads((root/'evidence.json').read_text())
    source=json.loads((Path(a.source)/'evidence.json').read_text());draft=json.loads(Path(a.draft).read_text())
    models={m['model_id']:m for m in data['models']};regions={r['region_id']:r for r in data['regions']}
    windows=fasta(root/'windows.fna');cds=fasta(root/'cds.fna');proteins=fasta(root/'proteins.faa')
    manifest={r['genome_id']:r for r in csv.DictReader(open(a.genomes),delimiter='\t')}
    checks=[];native={};start=time.perf_counter()
    def check(name,fn):
        t=time.perf_counter()
        try:details=fn();row=dict(name=name,status='passed',details=details)
        except Exception as error:row=dict(name=name,status='failed',error=type(error).__name__+': '+str(error))
        row['elapsed_seconds']=time.perf_counter()-t;checks.append(row);print(json.dumps(row),flush=True)
    def checkpoint():
        assert data['source_review_id']==source['review_id']==draft['review_id']
        assert data['source_evidence_sha256']==digest(Path(a.source)/'evidence.json')
        assert digest(root/'selection.input.json')==digest(a.draft)==data['selection']['decision_source_sha256']
        expected={l for pair in draft['pairs'] if pair.get('committed') is True for l in pair.get('selected',[])}
        assert {r['locus_id'] for r in regions.values()}==expected
        assert len(regions)==len(expected)==25 and len({r['genome_id'] for r in regions.values()})==13
        assert not data['selection'].get('explicit_genome_subset'),'Unexpected pilot restriction'
        lmap={l['locus_id']:l for l in source['loci']}
        for r in regions.values():
            l=lmap[r['locus_id']];assert (r['contig'],r['strand'],r['anchor_start'],r['anchor_end'])==(l['contig'],l['strand'],l['start'],l['end'])
            assert (r['window_start'],r['window_end'])==(max(0,l['start']-source['config']['window']),min(l['contig_length'],l['end']+source['config']['window']))
        assert set(windows)==set(regions) and len(models)==len(data['models'])
        return dict(regions=len(regions),genomes=13,omitted_unsaved_genomes=data['selection']['omitted_unsaved_genomes'])
    check('selection_and_native_window_bounds',checkpoint)
    def seals():
        files=json.loads((root/'seal.json').read_text())
        for name,h in files.items():
            path=(root/name).resolve();assert root.resolve() in path.parents and digest(path)==h,name
        json.dumps(data,allow_nan=False)
        for name in ['regions.decisions.tsv','models.decisions.tsv','neighbors.decisions.tsv']:
            for row in csv.DictReader(open(root/name),delimiter='\t'):assert not row.get('decision'),'Unexpected automatic decision'
        assert not list(root.glob('*msa*')) and not data['errors'],data['errors']
        return dict(sealed_files=len(files),no_automatic_approvals=True)
    check('seals_finite_values_and_pending_decisions',seals)
    def genomic_roundtrip():
        if a.reuse_native_audit:
            previous=json.loads(Path(a.reuse_native_audit).read_text());parent=data['provenance']['reidentification']
            parent_root=Path(parent['source']);original=json.loads((parent_root/'evidence.json').read_text())
            assert previous['status']=='passed' and previous['results']==str(parent_root) and previous['review_id']==parent['source_review_id']==original['review_id']
            assert parent['source_evidence_sha256']==digest(parent_root/'evidence.json')
            assert data['regions']==original['regions'] and data['provenance']['genomes']==original['provenance']['genomes']
            assert data['provenance']['prediction']==original['provenance']['prediction']
            assert [{k:v for k,v in m.items() if k!='identification'} for m in data['models']]==[{k:v for k,v in m.items() if k!='identification'} for m in original['models']]
            for name in ['windows.fna','cds.fna','proteins.faa','models.gff3','regions.tsv','annotations.tsv']:
                assert digest(root/name)==digest(parent_root/name),name
            for path in (root/'annotations').glob('*.gff3'):assert digest(path)==digest(parent_root/'annotations'/path.name)
            proof=next(c for c in previous['checks'] if c['name']=='independent_genomic_windows_CDS_and_protein_roundtrips')
            assert proof['status']=='passed'
            return dict(method='Reused independently verified native roundtrips after exact model/sequence/GFF/input/prediction identity checks',original_audit_sha256=digest(a.reuse_native_audit),original_proof=proof)
        proofs=[]
        for gid in sorted({r['genome_id'] for r in regions.values()}):
            entry=manifest[gid];path=Path(entry['fasta']);own=[r for r in regions.values() if r['genome_id']==gid]
            local=[m for m in models.values() if m['genome_id']==gid];needed={r['contig'] for r in own};seen=set()
            receipt=next(x for x in data['provenance']['genomes'] if x['genome_id']==gid)
            assert receipt['status']=='completed' and digest(path)==receipt['fasta_sha256']
            if entry.get('gff3'):assert digest(entry['gff3'])==receipt['annotation_sha256']
            with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as stream:
                for record in SeqIO.parse(stream,'fasta'):
                    if record.id not in needed:continue
                    assert record.id not in seen;seen.add(record.id);sequence=str(record.seq).upper()
                    for r in own:
                        if r['contig']==record.id:assert sequence[r['window_start']:r['window_end']]==windows[r['region_id']],r['region_id']
                    for m in local:
                        if m['contig']==record.id:native[m['model_id']]=native_cds(m,sequence)
                    if seen==needed:break
            assert seen==needed
            proofs.append(dict(genome_id=gid,windows=len(own),models=len(local),input_identity=receipt['identity_check']))
        assert set(native)==set(models)==set(cds)
        for mid,m in models.items():
            assert native[mid]==cds[mid],mid
            assert hashlib.sha256(cds[mid].encode()).hexdigest()==m['cds_sha256']
            aa=translation(cds[mid]).removesuffix('*')
            assert hashlib.sha256(aa.encode()).hexdigest()==m['protein_sha256'] and len(aa)==m['protein_length'],mid
            assert (mid in proteins)==bool(aa and '*' not in aa)
            if mid in proteins:assert proteins[mid]==aa
        return dict(genomes=proofs,model_roundtrips=len(models),usable_proteins=len(proteins),method='BioPython FASTA/strand handling/translation; original native genomic coordinates')
    check('independent_genomic_windows_CDS_and_protein_roundtrips',genomic_roundtrip)
    def gff_roundtrip():
        observed=defaultdict(list);stops=defaultdict(list);parents={}
        for line in (root/'models.gff3').read_text().splitlines():
            if line.startswith('#'):continue
            f=line.split('\t');attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
            if f[2]=='mRNA':parents[attrs['ID']]=attrs['Parent']
            if f[2] in {'CDS','stop_codon'}:
                mid=parents[attrs['Parent']];m=models[mid]
                assert f[6]==m['strand']
                assert tuple(unquote(x) for x in f[0].split('|'))==(m['genome_id'],m['contig'])
                if f[2]=='CDS':observed[mid].append([int(f[3])-1,int(f[4]),None if f[7]=='.' else int(f[7])])
                else:stops[mid].append([int(f[3])-1,int(f[4])])
        for mid,m in models.items():
            assert sorted(observed[mid])==sorted(m['cds_blocks'])
            assert sorted(stops[mid])==sorted(m['stop_blocks'])
        for row in csv.DictReader(open(root/'annotations.tsv'),delimiter='\t'):
            for line in (root/row['gff3']).read_text().splitlines():
                if line.startswith('#'):continue
                f=line.split('\t');attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
                if f[2]!='gene':continue
                model=models[attrs['ID']]
                assert row['genome_id']==model['genome_id'] and f[0]==model['contig']
        return dict(models=len(models),native_genome_gff_files=len(list((root/'annotations').glob('*.gff3'))))
    check('exported_GFF3_roundtrip',gff_roundtrip)
    def predictions():
        batches=data['provenance']['prediction']['batches']
        canonical={c['region_id']:c for c in data.get('prediction_intervals',[])}
        expected=set(canonical) if canonical else set(regions)
        assert len(batches)==len(expected)
        assert {rid for b in batches for rid in b['loci']}==expected
        canonical_dna=fasta(root/'prediction_windows.fna') if canonical else {}
        for b in batches:
            assert b['status']=='completed' and b['cpu_execution']['input_windows']==1
            assert not any('GPU' in s for s in b['cpu_execution']['visible_devices'])
            if canonical:
                cid=b['loci'][0];c=canonical[cid]
                assert len(canonical_dna[cid])==c['window_end']-c['window_start']
                for rid in c['region_ids']:
                    r=regions[rid]
                    assert canonical_dna[cid][r['window_start']-c['window_start']:r['window_end']-c['window_start']]==windows[rid]
                projected=[m for m in models.values() if m.get('prediction_interval_id')==cid]
                assert b['coding_models']==len({m['prediction_model_id'] for m in projected})
                raw=defaultdict(list)
                for line in (root/b['raw_gff3']).read_text().splitlines():
                    if line.startswith('#'):continue
                    f=line.split('\t')
                    if f[2]!='CDS':continue
                    attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
                    for parent in attrs['Parent'].split(','):
                        raw[parent].append([int(f[3])-1+c['window_start'],int(f[4])+c['window_start'],None if f[7]=='.' else int(f[7])])
                for m in projected:assert sorted(m['cds_blocks'])==sorted(raw[m['prediction_model_id']])
            else:
                assert b['coding_models']==sum(m['source']=='tiberius' and m['region_id'] in b['loci'] for m in models.values())
        return dict(batches=len(batches),coding_model_observations=sum(m['source']=='tiberius' for m in models.values()),
                    cpu_only=True,canonical_projection_checked=bool(canonical))
    check('all_prediction_workers_and_CPU_proofs',predictions)
    query_lengths={k:len(v) for k,v in proteins.items()};params=data['provenance']['parameters']
    def searches():
        db=Path(data['provenance']['protein_db_path']);info=data['provenance']['protein_db']
        for name,h in info['files'].items():assert digest(db/name)==h,name
        con=sqlite3.connect('file:'+str(db/'metadata.sqlite')+'?mode=ro',uri=True)
        refs={r[0]:r for r in con.execute('SELECT id,original_id,gene,length FROM proteins')};con.close()
        broad=raw_hits(root/'identification/blast.tsv',query_lengths,{k:v[3] for k,v in refs.items()})
        saved=json.loads((root/'identification_hits.json').read_text());assert {(h['query'],h['subject']) for h in saved}==set(broad)
        for h in saved:same_statistics(h,broad[h['query'],h['subject']])
        for m in models.values():
            if m['model_id'] not in proteins:continue
            passing=[(sid,h) for (qid,sid),h in broad.items() if qid==m['model_id'] and min(h['coverage_query'],h['coverage_subject'])>=params['identification_coverage'] and h['identity']>=params['identification_identity']]
            best=max((h['bitscore'] for _,h in passing),default=0)
            near={sid for sid,h in passing if h['bitscore']>=best*params['competing_ratio']}
            identification=m['identification'];assert identification['passing_hits']==len(passing)
            assert {h['subject'] for h in identification['near_best']}==near
            folded=identification.get('name_equivalence')=='symbol_case_only; accessions_exact'
            names={(('symbol',refs[sid][2].casefold()) if refs[sid][2] else ('accession',refs[sid][1])) if folded else (refs[sid][2] or refs[sid][1]) for sid in near}
            assert identification['status']==('no_passing_match' if not near else 'homology_assignment' if len(names)==1 else 'ambiguous_family_assignment')
            ordered=sorted(passing,key=lambda item:(-item[1]['bitscore'],item[1]['evalue'],item[0]))
            proposed=(refs[ordered[0][0]][2] or refs[ordered[0][0]][1]) if len(names)==1 else ''
            assert identification['proposed_name']==proposed
            if folded:assert identification['symbol_variants']==sorted({refs[sid][2] for sid in near if refs[sid][2]})
        targets={q['id']:len(q['sequence']) for q in source['queries']}
        assert {k:len(v) for k,v in fasta(root/'target_support/subjects.faa').items()}==targets
        support=raw_hits(root/'target_support/blast.tsv',query_lengths,targets)
        stored={(m['model_id'],h['subject']):h for m in models.values() for h in m['target_support']};assert set(stored)==set(support)
        for key,h in stored.items():same_statistics(h,support[key])
        return dict(reference_proteins=len(refs),database_pairs=len(broad),query_sequences=len(targets),target_support_pairs=len(support),database_release=info.get('download',{}).get('release'))
    check('raw_BLASTP_statistics_DB_naming_and_all_query_support',searches)
    def neighbors():
        local=raw_hits(root/'neighbor_support/blast.tsv',query_lengths,query_lengths);best={}
        for (a,b),h in local.items():
            ma,mb=models[a],models[b]
            if ma['genome_id']==mb['genome_id'] or min(h['coverage_query'],h['coverage_subject'])<params['identification_coverage'] or h['identity']<params['identification_identity']:continue
            key=ma['gene_group'],mb['gene_group']
            if key not in best or h['bitscore']>best[key]['bitscore']:best[key]=dict(h,a=a,b=b)
        near=defaultdict(set)
        for (ga,gb),h in best.items():
            group=ga,models[h['b']]['genome_id'];peers=[v['bitscore'] for (x,y),v in best.items() if x==ga and models[v['b']]['genome_id']==group[1]]
            if h['bitscore']>=max(peers)*params['competing_ratio']:near[group].add(gb)
        expected={(ga,gb) for (ga,gid),values in near.items() for gb in values}
        assert {(p['gene_a'],p['gene_b']) for p in data['neighbor_pairs']}==expected
        focal={m['gene_group'] for m in models.values() if m['overlaps_anchor']}
        for p in data['neighbor_pairs']:
            ma,mb=models[p['model_a']],models[p['model_b']];assert ma['genome_id']!=mb['genome_id']
            assert (ma['gene_group'],mb['gene_group'])==(p['gene_a'],p['gene_b'])
            # Equal-score isoforms may have different lengths/coverage. Check
            # the actual selected model pair, without imposing a tie order.
            same_statistics(p,local[p['model_a'],p['model_b']])
            assert p['bitscore']==best[p['gene_a'],p['gene_b']]['bitscore']
            forward=near[p['gene_a'],mb['genome_id']];reverse=near[p['gene_b'],ma['genome_id']]
            unique=len(forward)==len(reverse)==1 and p['gene_a'] in reverse
            assert p['status']==('unique_reciprocal_proposal' if unique else 'ambiguous_proposal')
            assert p['kind']==('target_similarity' if p['gene_a'] in focal or p['gene_b'] in focal else 'neighbor_proposal')
            assert p['orthology_confirmed'] is False
        return dict(raw_model_pairs=len(local),directed_gene_proposals=len(expected),unordered_gene_proposals=len({tuple(sorted(v)) for v in expected}))
    check('isoform_collapse_neighbor_competition_and_target_exclusion',neighbors)
    if a.pilot:
        def replay():
            pilot=json.loads((Path(a.pilot)/'evidence.json').read_text());count=0
            for m in pilot['models']:
                current=models[m['model_id']]
                for key in ['start','end','strand','cds_blocks','stop_blocks','cds_sha256','protein_sha256','source']:assert current[key]==m[key],(m['model_id'],key)
                count+=1
            return dict(unchanged_pilot_model_roundtrips=count)
        check('independent_window_replay_against_original_pilot',replay)
    rows=[]
    for r in regions.values():
        own=[m for m in models.values() if m['region_id']==r['region_id']];anchor=[m for m in own if m['overlaps_anchor']]
        best=sorted([(m,h) for m in anchor if m['source']=='tiberius' for h in m['target_support']],key=lambda v:-v[1]['bitscore'])
        rows.append(dict(r,annotation_models=sum(m['source']=='annotation' for m in own),predicted_models=sum(m['source']=='tiberius' for m in own),
            annotation_anchor_models=sum(m['source']=='annotation' for m in anchor),predicted_anchor_models=sum(m['source']=='tiberius' for m in anchor),
            models_with_issues=sum(bool(m['issues']) for m in own),edge_models=sum(m['window_edge'] for m in own),
            best_anchor_prediction_support=dict(model=best[0][0]['model_id'],**best[0][1]) if best else None,
            annotation_names=sorted({m['annotation_name'] for m in anchor if m['annotation_name']})))
    output=dict(scope='All 25 saved candidate loci; 13 species; no biological model approval',status='passed' if all(c['status']=='passed' for c in checks) else 'failed',
        results=str(root),review_id=data['review_id'],checks=checks,elapsed_seconds=time.perf_counter()-start,regions=rows,
        totals=dict(regions=len(regions),genomes=len({r['genome_id'] for r in regions.values()}),model_alternatives=len(models),model_sources=dict(Counter(m['source'] for m in models.values())),
            coding_status=dict(Counter(m['coding_status'] for m in models.values())),identification_status=dict(Counter(m['identification']['status'] for m in models.values()))),
        limitations=['Structural and software verification does not establish expression, function, orthology or biological completeness.','The broad reviewed vertebrate database is not an exhaustive proteome.','A missing ab initio prediction does not establish gene absence.'])
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print('AUDIT '+output['status'].upper()+' '+str(out),flush=True)
    return 0 if output['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
