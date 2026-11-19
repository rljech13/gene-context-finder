"""One durable host-local queue for stage-2 commands, including live adoption.

The dispatcher alone admits work. Scientific commands and checkpoints stay pinned.
Linux /proc identities prevent PID reuse from resuming an unrelated process.
"""
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import html
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

TERMINAL={'awaiting_stop2_review','needs_context_review','completed_with_errors','failed'}
DEFAULTS=dict(max_active=2,max_threads=64,max_job_threads=32,
              memory_max_bytes=250*1024**3,memory_high_bytes=200*1024**3,
              admission_memory_bytes=192*1024**3,reserve_bytes=32*1024**3,poll_seconds=2)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()


def atomic(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.pending.'+str(os.getpid()))
    with tmp.open('w') as f:
        json.dump(value,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)


def process(pid):
    try:
        p=Path('/proc')/str(pid);fields=(p/'stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=int(pid),ppid=int(fields[1]),state=fields[0],start_ticks=int(fields[19]),
                    wait_status=int(fields[49]) if len(fields)>49 else None,
                    argv=[v.decode() for v in (p/'cmdline').read_bytes().split(b'\0') if v])
    except (OSError,ValueError,ProcessLookupError):return None


def matches(ref):
    p=process(ref['pid'])
    return p if p and p['start_ticks']==ref['start_ticks'] else None


def descendants(pid):
    values={}
    for entry in Path('/proc').iterdir():
        if entry.name.isdigit():
            p=process(int(entry.name))
            if p:values[p['pid']]=p
    own={pid}
    while True:
        expanded=own|{n for n,p in values.items() if p['ppid'] in own}
        if own==expanded:break
        own=expanded
    return [values[n] for n in sorted(own) if n in values]


def pause_tree(ref):
    if not matches(ref):return []
    os.kill(ref['pid'],signal.SIGSTOP)
    for _ in range(3):
        own=descendants(ref['pid'])
        for p in own:
            if p['state']!='Z':
                try:os.kill(p['pid'],signal.SIGSTOP)
                except ProcessLookupError:pass
    return own


def resume_tree(ref):
    if not matches(ref):raise RuntimeError('Process identity changed before resume')
    for p in descendants(ref['pid']):
        if p['state']!='Z':
            try:os.kill(p['pid'],signal.SIGCONT)
            except ProcessLookupError:pass


def identity(spec):
    """Ignore output placement/CPU width; bind immutable inputs and science args."""
    args=spec['command'];canonical=[];i=0
    while i<len(args):
        value=args[i]
        if value in {'--out','--threads'}:i+=2;continue
        if value in {'--draft','--decisions','--genomes','--prediction-models'}:
            canonical.extend([value,digest(args[i+1])]);i+=2;continue
        if value=='--results':
            canonical.extend([value,digest(Path(args[i+1])/'evidence.json')]);i+=2;continue
        canonical.append(value);i+=1
    return hashlib.sha256(json.dumps([spec['source'],canonical],sort_keys=True).encode()).hexdigest()


def admit(policy,jobs,threads,nonreclaimable):
    active=[j for j in jobs if j['status'] in {'running','launching'}]
    return (len(active)<policy['max_active'] and
            sum(j['threads'] for j in active)+threads<=policy['max_threads'] and
            nonreclaimable+policy['reserve_bytes']+policy.get('active_growth_reserve_bytes',policy['reserve_bytes'])*len(active)<=policy['admission_memory_bytes'])


def memory_state(cgroup):
    cgroup=Path(cgroup)
    stats={k:int(v) for k,v in (line.split() for line in (cgroup/'memory.stat').read_text().splitlines())}
    # file cache can be reclaimed; include dirty/writeback pages conservatively.
    nonreclaimable=sum(stats.get(k,0) for k in ['anon','slab_unreclaimable','kernel_stack','pagetables','file_dirty','file_writeback','shmem'])
    return dict(current_bytes=int((cgroup/'memory.current').read_text()),nonreclaimable_bytes=nonreclaimable,
                file_cache_bytes=stats.get('file',0),hard_limit=(cgroup/'memory.max').read_text().strip(),
                swap_limit=(cgroup/'memory.swap.max').read_text().strip(),soft_limit=(cgroup/'memory.high').read_text().strip(),
                events={k:int(v) for k,v in (line.split() for line in (cgroup/'memory.events').read_text().splitlines())})


class Queue:
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)

    @contextmanager
    def locked(self):
        with (self.root/'state.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            state=json.loads((self.root/'state.json').read_text())
            yield state
            atomic(self.root/'state.json',state)

    def event(self,kind,**fields):
        with (self.root/'events.jsonl').open('a') as f:
            f.write(json.dumps(dict(at=time.time(),kind=kind,**fields),ensure_ascii=False)+'\n');f.flush();os.fsync(f.fileno())

    def add(self,spec):
        argv=spec['command']
        if argv[1:4]!=['-m','genecontext','stage2']:raise ValueError('Only a pinned stage2 command is accepted')
        threads=int(argv[argv.index('--threads')+1]);key=identity(spec)
        with self.locked() as state:
            if not 1<=threads<=state['policy']['max_job_threads']:raise ValueError('Per-job thread budget exceeded')
            for j in state['jobs']:
                if j['identity']==key:return j['id'],False
            j=dict(spec,id='Q'+str(len(state['jobs'])+1).zfill(4),identity=key,threads=threads,
                   enqueued=time.time(),status='queued',process=None,adopted=False)
            out=Path(argv[argv.index('--out')+1])
            if out.exists():raise ValueError('Output already exists; adopt explicitly instead of rerunning')
            j['output']=str(out)
            state['jobs'].append(j);self.event('enqueued',id=j['id'],panel=j['panel'])
            return j['id'],True

    def render(self,state):
        delivery=self.root/'delivery';delivery.mkdir(exist_ok=True)
        for name in ['state.json','events.jsonl']:
            link=delivery/name
            if not link.exists():link.symlink_to(self.root/name)
        esc=lambda v:html.escape(str(v))
        labels={'finishing':'Проверка готового результата','launching':'Запускается','queued':'В очереди','paused':'В очереди — прогресс сохранён','running':'Выполняется',
                'awaiting_stop2_review':'Готово к STOP 2','needs_context_review':'STOP 2: проверить окружение',
                'completed_with_errors':'STOP 2: есть ошибки','failed':'Ошибка — требуется разбор'}
        rows=[];pending=0
        for j in state['jobs']:
            if j['status'] in {'queued','paused'}:pending+=1;place=str(pending)
            else:place='—'
            links='Ожидается'
            if j.get('published') and j.get('archive_link'):
                links=f"<a href='results/{j['id']}/report.html'>Отчёт</a> · <a href='exports/{esc(j['archive_link'])}'>ZIP</a>"
            rows.append(f"<tr><td>{j['id']}</td><td>{esc(j['panel'])}</td><td>{esc(j['name'])}</td><td>{j['loci']}</td><td>{place}</td><td>{labels[j['status']]}</td><td>{links}</td></tr>")
        m=state.get('memory',{});anon=m.get('nonreclaimable_bytes',0)/1024**3
        policy=state.get('policy',DEFAULTS)
        limits=f"Не более {policy['max_active']} анализов одновременно, суммарно до {policy['max_threads']} потоков."
        growth=policy.get('active_growth_reserve_bytes',policy['reserve_bytes'])/1024**3
        reserve=policy['reserve_bytes']/1024**3
        pause=policy.get('pause_memory_bytes',policy['memory_high_bytes'])/1024**3
        memory_note=f'Для нового запуска требуется резерв {reserve:g} GiB плюс по {growth:g} GiB на рост каждой работающей задачи. При невытесняемой памяти {pause:g} GiB работа приостанавливается.'
        page="""<!doctype html><html lang='ru'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><meta http-equiv='refresh' content='15'><title>Общая очередь этапа 2</title><style>body{font:16px/1.5 system-ui;max-width:1250px;margin:30px auto;padding:0 20px;color:#203044}table{border-collapse:collapse;width:100%}td,th{padding:10px;text-align:left;border-bottom:1px solid #cbd5df}.note{background:#eaf3fb;padding:16px}a{color:#245d9c}</style><h1>Общая очередь этапа 2</h1><p class='note'>Единый порядок для MYC, KLF4, KLF5 и следующих добавлений. QUEUE_CPU_LIMITS Общий жёсткий предел 250 GiB RAM, swap выключен; мягкий порог памяти 200 GiB. QUEUE_MEMORY_POLICY Полные научные окна и ручные группы сохраняются.</p><p>Приостановленные расчёты возобновляются с сохранённого состояния процесса. Готовые результаты не пересчитываются. После завершения — ручной STOP 2, без автоматического одобрения моделей и этапа 3.</p>"""+f"<p>Работают: {sum(j['status']=='running' for j in state['jobs'])}. В очереди: {pending}. Невытесняемая память общей группы: {anon:.1f} GiB; файловый кеш учитывается отдельно. <a href='state.json'>Статус JSON</a></p><table><tr><th>ID</th><th>Панель</th><th>Группа / вход</th><th>Локусы</th><th>Место</th><th>Статус</th><th>Результат</th></tr>"+''.join(rows)+'</table></html>'
        tmp=delivery/'index.pending.html';tmp.write_text(page.replace('QUEUE_CPU_LIMITS',limits).replace('QUEUE_MEMORY_POLICY',memory_note));os.replace(tmp,delivery/'index.html')

    def publish(self,j,code):
        """Audit/publish one completed checkpoint without repeating inference."""
        from collections import Counter,defaultdict
        import csv,zipfile
        target=Path(j['command'][j['command'].index('--out')+1])
        if not (target/'seal.json').exists():
            raise RuntimeError('Command ended without a sealed checkpoint; returncode='+str(code))
        evidence=json.loads((target/'evidence.json').read_text());seal=json.loads((target/'seal.json').read_text())
        assert all(digest(target/name)==h for name,h in seal.items())
        assert {r['locus_id'] for r in evidence['regions']}==set(j['locus_ids'])
        assert evidence['source_review_id']==j['source_review_id']
        assert (evidence['selection'].get('analysis_group') or {}).get('group_id')==j.get('group_id')
        for name in ['models.decisions.tsv','regions.decisions.tsv','neighbors.decisions.tsv']:
            assert all(not row.get('decision') for row in csv.DictReader((target/name).open(),delimiter='\t'))
        for name,h in j['source_files'].items():assert digest(Path(j['source'])/name)==h
        checks=[]
        for r in evidence['regions']:
            genes=defaultdict(list)
            for m in evidence['models']:
                if m['region_id']==r['region_id'] and m['source']=='tiberius':genes[(m.get('prediction_interval_id'),m['original_gene'])].append(m)
            sides=dict(left=0,right=0)
            for models in genes.values():
                if any(m['window_edge'] for m in models):continue
                a=min(m['start'] for m in models);b=max(m['end'] for m in models)
                if b<=r['anchor_start']:sides['left']+=1
                elif a>=r['anchor_end']:sides['right']+=1
            checks.append(dict(locus_id=r['locus_id'],predicted_distinct_genes=sides,minimum_five_each_side=min(sides.values())>=5))
        warnings=sum(not c['minimum_five_each_side'] for c in checks)
        status='completed_with_errors' if code or evidence['errors'] else 'needs_context_review' if warnings else 'awaiting_stop2_review'
        proof=dict(status=status,queue_id=j['id'],seal_verified=True,selected_locus_ids_exact=True,decisions_pending=True,
                   returncode=code,predicted_context=checks,context_warnings=warnings,errors=evidence['errors'],
                   source_evidence_sha256=digest(target/'evidence.json'),report_sha256=digest(target/'report.html'))
        atomic(target/'queue.audit.json',proof)
        archive=Path(j['legacy_batch']).parent/'exports'/(j.get('archive') or j['job_id']+'_STOP2.zip') if j.get('legacy_batch') else self.root/'exports'/(j['id']+'_STOP2.zip')
        archive.parent.mkdir(exist_ok=True)
        if not archive.exists():
            names=['report.html','evidence.json','window_plan.json','models.gff3','annotations.tsv','models.tsv','regions.tsv','gene_identities.tsv','cds.fna','proteins.faa','models.decisions.tsv','regions.decisions.tsv','neighbors.decisions.tsv','identification_hits.json','queue.audit.json','selection.snapshot.tsv','selection.input.json']
            files={n:target/n for n in names if (target/n).is_file()};files.update({p.relative_to(target).as_posix():p for p in (target/'annotations').glob('*.gff3')})
            manifest=dict(schema='gcf_queue_stop2_delivery_v1',queue_id=j['id'],files={n:dict(sha256=digest(p),bytes=p.stat().st_size) for n,p in files.items()},human_model_approval=False)
            tmp=archive.with_suffix('.zip.pending')
            with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
                for name,p in files.items():z.write(p,j['job_id']+'/'+name)
                z.writestr(j['job_id']+'/CHECKSUMS.json',json.dumps(manifest,indent=2)+'\n')
            os.replace(tmp,archive)
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None
            manifest=json.loads(z.read(j['job_id']+'/CHECKSUMS.json'))
            for name,v in manifest['files'].items():assert hashlib.sha256(z.read(j['job_id']+'/'+name)).hexdigest()==v['sha256']
        delivery=self.root/'delivery';(delivery/'results').mkdir(exist_ok=True);(delivery/'exports').mkdir(exist_ok=True)
        link=delivery/'results'/j['id']
        if not link.exists():link.symlink_to(target,target_is_directory=True)
        archive_link=j['id']+'_'+archive.name;link=delivery/'exports'/archive_link
        if not link.exists():link.symlink_to(archive)
        j.update(status=status,published=True,review_id=evidence['review_id'],archive=str(archive),archive_link=archive_link,
                 archive_sha256=digest(archive),report_sha256=proof['report_sha256'],evidence_sha256=proof['source_evidence_sha256'],
                 models=dict(Counter(m['source'] for m in evidence['models'])),context_warnings=warnings,errors=evidence['errors'])

    def sync_legacy(self,state):
        for path in {j.get('legacy_batch') for j in state['jobs'] if j.get('legacy_batch')}:
            p=Path(path);batch=json.loads(p.read_text())
            for j in state['jobs']:
                if j.get('legacy_batch')!=path:continue
                g=next(g for g in batch['groups'] if g['job_id']==j['job_id'])
                g.update(status=j['status'],queue_id=j['id'],shared_queue=True)
                for key in ['published','review_id','archive_sha256','report_sha256','evidence_sha256','models','errors','returncode']:
                    if key in j:g[key]=j[key]
                if j.get('archive'):g['archive']=Path(j['archive']).name
            own=[j for j in state['jobs'] if j.get('legacy_batch')==path]
            batch['status']='completed_for_stop2_review' if all(j['status'] in TERMINAL and j.get('published') for j in own) else 'completed_with_errors' if all(j['status'] in TERMINAL for j in own) else 'managed_by_shared_queue'
            batch['shared_queue_url']='/stage2_queue/';atomic(p,batch)

    def retire_parents(self,state):
        for parent in state.get('legacy_parents',[]):
            if parent.get('retired'):continue
            if any(j['adopted'] and j.get('legacy_batch')==parent['batch_path'] and j['status'] not in TERMINAL for j in state['jobs']):continue
            p=matches(parent)
            if p:
                own=descendants(p['pid'])
                if any(v['pid']!=p['pid'] and v['state']!='Z' for v in own):continue
                os.kill(p['pid'],signal.SIGKILL)
            parent['retired']=time.time();self.event('legacy_dispatcher_retired',pid=parent['pid'],reason='All adopted scientific children completed; no pending old commands were allowed to start')

    def step(self):
        with self.locked() as state:
            policy=state['policy'];memory=memory_state(state['cgroup'])
            if memory['hard_limit']!=str(policy['memory_max_bytes']) or memory['swap_limit']!='0' or memory['soft_limit']!=str(policy['memory_high_bytes']):raise RuntimeError('Hard shared memory/swap policy changed; dispatch refused')
            state['memory']=memory
            oom=memory['events'].get('oom_kill',0)>state.get('oom_kill_baseline',0)
            pressure=memory['nonreclaimable_bytes']>=policy.get('pause_memory_bytes',policy['memory_high_bytes'])
            state['memory_hold']=bool(oom or pressure)
            for j in state['jobs']:
                if j['status'] in {'running','paused'} and j.get('process'):
                    actual=matches(j['process'])
                    if actual and actual['state'] not in {'Z','X'}:
                        j['status']='paused' if actual['state'] in {'T','t'} else 'running'
                    if state['memory_hold'] and j['status']=='running':
                        pause_tree(j['process']);j['status']='paused';self.event('memory_pause',id=j['id'],oom_event=oom)
            for j in state['jobs']:
                if j['status'] not in {'running','paused','finishing','launching'}:continue
                if j['status']=='launching':
                    receipt=self.root/'jobs'/j['id']/'process.json'
                    if receipt.exists():
                        j['process']=json.loads(receipt.read_text());j['status']='running'
                    elif time.time()-j['launch_requested']<30:continue
                    else:
                        j.update(status='failed',error='Launch did not produce a process receipt; no automatic duplicate/retry')
                        self.event('launch_failed',id=j['id']);continue
                p=matches(j['process']) if j.get('process') else None
                exit_path=self.root/'jobs'/j['id']/'exit.json'
                code=None
                if exit_path.exists():code=json.loads(exit_path.read_text())['returncode']
                elif p and p['state']=='Z':code=os.waitstatus_to_exitcode(p['wait_status'])
                elif not p:code=0 if (Path(j['output'])/'seal.json').exists() else -1
                if code is not None:
                    j.update(status='finishing',returncode=code,finished=time.time())
                    try:self.publish(j,code)
                    except Exception as e:j.update(status='failed',error=str(e))
                    self.event('finished',id=j['id'],status=j['status'],returncode=code)
            # A worker restart does not start another copy of any live process.
            active=[j for j in state['jobs'] if j['status']=='running']
            for j in reversed(active):
                if len([x for x in state['jobs'] if x['status']=='running'])<=policy['max_active'] and sum(x['threads'] for x in state['jobs'] if x['status']=='running')<=policy['max_threads']:break
                pause_tree(j['process']);j['status']='paused';self.event('paused',id=j['id'])
            for j in state['jobs']:
                if j['status'] not in {'queued','paused'}:continue
                memory=memory_state(state['cgroup']);state['memory']=memory
                if state['memory_hold'] or not admit(policy,state['jobs'],j['threads'],memory['nonreclaimable_bytes']):break
                if time.time()-state.get('last_admission',0)<policy.get('admission_interval_seconds',0):break
                state['last_admission']=time.time()
                if j['status']=='paused':
                    resume_tree(j['process']);self.event('resumed',id=j['id'])
                else:
                    for name,h in j['source_files'].items():
                        if digest(Path(j['source'])/name)!=h:raise RuntimeError('Pinned scientific source changed')
                    jobdir=self.root/'jobs'/j['id'];jobdir.mkdir(parents=True,exist_ok=True);atomic(jobdir/'spec.json',j)
                    j['status']='launching';j['launch_requested']=time.time();atomic(self.root/'state.json',state)
                    # Each scientific job has its own bounded systemd service.
                    # Restarting the dispatcher cannot terminate these children.
                    argv=[sys.executable,state['submit_script'],'--name','stage2-'+j['id'].lower(),'--threads',str(j['threads']),'--',sys.executable,str(Path(__file__).resolve()),'child','--spec',str(jobdir/'spec.json'),'--exit',str(jobdir/'exit.json')]
                    result=subprocess.run(argv,capture_output=True,text=True)
                    if result.returncode:raise RuntimeError('Job submission failed: '+result.stderr[-1000:])
                    submission=json.loads(result.stdout);atomic(jobdir/'submission.json',submission)
                    j['service']=submission['unit'];self.event('submitted',id=j['id'],service=j['service'])
                    continue
                j['status']='running'
            self.retire_parents(state);self.sync_legacy(state);state['updated']=time.time();self.render(state)

    def run(self):
        with (self.root/'worker.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            while True:
                try:self.step()
                except Exception as e:
                    self.event('dispatcher_error',error=repr(e));print(repr(e),file=sys.stderr,flush=True)
                time.sleep(DEFAULTS['poll_seconds'])


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    for name in ['worker','status','enqueue']:
        c=sub.add_parser(name);c.add_argument('--root',required=True)
        if name=='enqueue':c.add_argument('--spec',required=True)
    c=sub.add_parser('child');c.add_argument('--spec',required=True);c.add_argument('--exit',required=True)
    args=p.parse_args()
    if args.action=='child':
        spec_path=Path(args.spec)
        with (spec_path.parent/'execution.lock').open('a') as execution_lock:
            fcntl.flock(execution_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            if Path(args.exit).exists():return json.loads(Path(args.exit).read_text())['returncode']
            ref=process(os.getpid());atomic(spec_path.parent/'process.json',dict(pid=ref['pid'],start_ticks=ref['start_ticks']))
            spec=json.loads(spec_path.read_text());env=dict(os.environ,PYTHONPATH=spec['source'],GCF_THREADS=str(spec['threads']),OMP_NUM_THREADS=str(spec['threads']),MKL_NUM_THREADS=str(spec['threads']),PYTHONUNBUFFERED='1')
            result=subprocess.run(spec['command'],cwd=spec['source'],env=env)
            atomic(args.exit,dict(returncode=result.returncode,finished=time.time(),command=spec['command']))
            return result.returncode
    queue=Queue(args.root)
    if args.action=='worker':queue.run()
    elif args.action=='enqueue':
        identifier,created=queue.add(json.loads(Path(args.spec).read_text()));print(json.dumps(dict(id=identifier,created=created)))
    else:print((queue.root/'state.json').read_text())
    return 0


if __name__=='__main__':raise SystemExit(main())
