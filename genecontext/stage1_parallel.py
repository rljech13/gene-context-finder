"""Independent genome workers under the caller's shared resource boundary."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import os
import shutil
import subprocess
import sys
from uuid import uuid4

from .pipeline import load_genomes
from .sequence import sha256, write_tsv
from .stage1 import dump, new_output, now, verify, write_outputs


OPTIONS = ['gene','query','query_type','query_metadata','threads','search_evalue',
           'seed_min_length','seed_min_identity','max_intron','chain_max_overlap',
           'max_locus_span','locus_shared_fraction','window','refine','refine_padding',
           'refine_max_models','min_query_coverage','min_identity','max_evalue']


def merge_parts(parts, output, expected_genomes):
    """Verify worker checkpoints and produce one ordinary stage 1 checkpoint."""
    data=None;seen=set();receipts=[];roots=[]
    for part in parts:
        root,current=verify(part);roots.append(root)
        if data is None:
            data=dict(current,loci=[],genomes=[],commands=[])
        if (current['queries']!=data['queries'] or current['config']!=data['config'] or
                current['filters']!=data['filters'] or current['query_provenance']!=data['query_provenance']):
            raise ValueError('Cannot merge workers with different queries/settings')
        ids={g['genome_id'] for g in current['genomes']}
        if seen & ids:raise ValueError('Duplicate genome IDs in worker checkpoints')
        seen.update(ids);data['genomes'].extend(current['genomes']);data['loci'].extend(current['loci'])
        data['commands'].extend(current.get('commands',[]))
        receipts.append(dict(path=str(root),review_id=current['review_id'],evidence_sha256=sha256(root/'evidence.json')))
    if data is None:raise ValueError('No worker checkpoints to merge')
    if seen!=set(expected_genomes):raise ValueError('Missing/unexpected genome worker results')
    order={gid:i for i,gid in enumerate(expected_genomes)}
    data['genomes'].sort(key=lambda g:order[g['genome_id']])
    data['loci'].sort(key=lambda l:(order[l['genome_id']],l['contig'],l['start'],l['locus_id']))
    data['review_id']=str(uuid4());data['source_part_reviews']=receipts
    root=new_output(output);(root/'raw').mkdir()
    shutil.copyfile(roots[0]/'queries.fasta',root/'queries.fasta')
    with (root/'loci.fna').open('wb') as out:
        for part in roots:
            if (part/'loci.fna').exists():
                with (part/'loci.fna').open('rb') as inp:shutil.copyfileobj(inp,out)
            for path in (part/'raw').iterdir():
                target=root/'raw'/path.name
                if target.exists():raise ValueError('Worker raw evidence names collide: '+path.name)
                if path.is_dir():shutil.copytree(path,target)
                else:shutil.copyfile(path,target)
    write_outputs(root,data)
    return data


def search(args):
    # Only processes are parallel: every worker runs the same stage1 algorithm.
    genomes=load_genomes(args)
    root=new_output(args.out);work=root/'workers';work.mkdir()
    rows={g['genome_id']:dict(g) for g in genomes}
    statuses={gid:dict(genome_id=gid,status='queued') for gid in rows}
    options=[]
    for key in OPTIONS:
        value=getattr(args,key)
        if value is not None:
            if key in {'query','query_metadata'}:value=str(Path(value).resolve())
            options+=['--'+key.replace('_','-'),str(value)]
    def worker(row):
        gid=row['genome_id'];folder=work/gid;folder.mkdir()
        manifest=folder/'genome.tsv';write_tsv(manifest,[row],list(row))
        output=folder/'result'
        argv=[sys.executable,'-m','genecontext','stage1',*options,'--jobs','1',
              '--genomes',str(manifest),'--out',str(output)]
        env=dict(os.environ,OMP_NUM_THREADS=str(args.threads))
        with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
            result=subprocess.run(argv,stdout=out,stderr=err,env=env)
        dump(folder/'completed.json',dict(argv=argv,returncode=result.returncode,finished=now()))
        return output,result.returncode
    futures={};completed={};failures=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for row in genomes:futures[pool.submit(worker,row)]=row['genome_id']
        for future in as_completed(futures):
            gid=futures[future]
            try:
                output,code=future.result();verify(output);completed[gid]=output
                statuses[gid].update(status='completed' if code==0 else 'completed_with_errors',returncode=code)
                if code:failures.append(gid)
            except (OSError,ValueError,RuntimeError) as error:
                statuses[gid].update(status='worker_failed',error=str(error));failures.append(gid)
            for pending,name in futures.items():
                if not pending.done():statuses[name]['status']='running' if pending.running() else 'queued'
            dump(root/'progress.json',dict(updated=now(),jobs=args.jobs,threads_per_genome=args.threads,
                                           genomes=list(statuses.values())))
            print(f"[stage 1] {len(completed)}/{len(genomes)} genomes saved: {gid}",file=sys.stderr,flush=True)
    if len(completed)!=len(genomes):
        dump(root/'batch_failed.json',dict(failures=failures,genomes=list(statuses.values())))
        raise RuntimeError('Some workers did not produce valid checkpoints; see '+str(root/'progress.json'))
    # Aggregate into the output root so the normal report/selection paths work.
    merged=root/'merged'
    data=merge_parts([completed[g['genome_id']] for g in genomes],merged,[g['genome_id'] for g in genomes])
    for path in merged.iterdir():shutil.move(str(path),root/path.name)
    merged.rmdir()
    dump(root/'batch.json',dict(jobs=args.jobs,threads_per_genome=args.threads,
                               source_part_reviews=data['source_part_reviews'],failures=failures))
    print('STOP — open '+str(root/'report.html'),file=sys.stderr)
    return 1 if failures else 0
