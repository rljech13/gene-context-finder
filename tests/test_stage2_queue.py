from copy import deepcopy
import json,os,signal,subprocess,sys,time
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from genecontext.stage2_queue import Queue,DEFAULTS,admit,atomic,identity,matches,pause_tree,process,resume_tree

class SharedQueueTests(unittest.TestCase):
    def test_common_limits_apply_across_panels_and_paused_jobs_do_not_take_cpu_slot(self):
        jobs=[dict(panel='MYC',status='running',threads=32),dict(panel='KLF4',status='running',threads=32),dict(panel='KLF5',status='paused',threads=16)]
        self.assertFalse(admit(DEFAULTS,jobs,16,20*1024**3))
        jobs[0]['status']='awaiting_stop2_review'
        self.assertTrue(admit(DEFAULTS,jobs,16,20*1024**3))
        self.assertFalse(admit(DEFAULTS,[dict(status='running',threads=48)],32,0))

    def test_four_job_policy_keeps_memory_reserve_and_dashboard_consistent(self):
        policy=dict(DEFAULTS,max_active=4,max_threads=128)
        jobs=[dict(status='running',threads=32) for _ in range(3)]
        self.assertTrue(admit(policy,jobs,32,20*1024**3))
        self.assertFalse(admit(policy,jobs,32,65*1024**3))
        self.assertFalse(admit(policy,jobs+[dict(status='running',threads=32)],32,0))
        with TemporaryDirectory() as t:
            q=Queue(t);q.render(dict(policy=policy,jobs=[]))
            page=(q.root/'delivery/index.html').read_text()
            self.assertIn('Не более 4 анализов одновременно, суммарно до 128 потоков.',page)
            self.assertNotIn('QUEUE_CPU_LIMITS',page)

    def test_eight_jobs_require_startup_and_growth_reserves(self):
        policy=dict(DEFAULTS,max_active=8,max_threads=256,active_growth_reserve_bytes=16*1024**3,pause_memory_bytes=128*1024**3,admission_interval_seconds=30)
        jobs=[dict(status='running',threads=32) for _ in range(7)]
        self.assertTrue(admit(policy,jobs,32,20*1024**3))
        self.assertFalse(admit(policy,jobs,32,49*1024**3))
        self.assertFalse(admit(policy,jobs+[dict(status='running',threads=32)],16,0))
        with TemporaryDirectory() as t:
            q=Queue(t);q.render(dict(policy=policy,jobs=[]))
            page=(q.root/'delivery/index.html').read_text()
            self.assertIn('Не более 8 анализов одновременно, суммарно до 256 потоков.',page)
            self.assertIn('резерв 32 GiB плюс по 16 GiB',page)
            self.assertIn('невытесняемой памяти 128 GiB',page)
            self.assertNotIn('QUEUE_MEMORY_POLICY',page)

    def test_memory_reserve_blocks_new_work_well_before_hard_limit(self):
        self.assertTrue(admit(DEFAULTS,[],32,159*1024**3))
        self.assertFalse(admit(DEFAULTS,[],32,161*1024**3))
        self.assertLess(DEFAULTS['admission_memory_bytes'],DEFAULTS['memory_high_bytes'])
        self.assertLess(DEFAULTS['memory_high_bytes'],DEFAULTS['memory_max_bytes'])

    def test_enqueue_is_idempotent_for_relocated_outputs_and_preserves_fifo(self):
        with TemporaryDirectory() as t:
            root=Path(t);queue=Queue(root/'queue');atomic(queue.root/'state.json',dict(policy=DEFAULTS,jobs=[]))
            source=root/'source';source.mkdir();(source/'evidence.json').write_text('{}')
            draft=root/'draft.json';draft.write_text('{"selection":1}')
            spec=dict(source='pinned-source',panel='KLF5',command=[sys.executable,'-m','genecontext','stage2','--results',str(source),'--draft',str(draft),'--out',str(root/'a'),'--threads','16'])
            self.assertEqual(queue.add(spec),('Q0001',True))
            other=deepcopy(spec);other['command'][other['command'].index('--out')+1]=str(root/'b')
            self.assertEqual(queue.add(other),('Q0001',False))
            draft2=root/'different.json';draft2.write_text('{"selection":2}');other['command'][other['command'].index('--draft')+1]=str(draft2)
            self.assertEqual(queue.add(other),('Q0002',True))
            saved=json.loads((queue.root/'state.json').read_text());self.assertEqual([j['id'] for j in saved['jobs']],['Q0001','Q0002'])

    def test_refuses_unmanaged_command_and_excess_threads(self):
        with TemporaryDirectory() as t:
            q=Queue(t);atomic(q.root/'state.json',dict(policy=DEFAULTS,jobs=[]))
            with self.assertRaisesRegex(ValueError,'Only a pinned'):q.add(dict(command=['echo','wrong']))
            with self.assertRaisesRegex(ValueError,'thread budget'):q.add(dict(source='source',command=[sys.executable,'-m','genecontext','stage2','--out',str(Path(t)/'out'),'--threads','65']))

    def test_migration_renders_completed_legacy_result_before_new_links_exist(self):
        with TemporaryDirectory() as t:
            q=Queue(t)
            q.render(dict(jobs=[dict(id='Q0001',panel='MYC',name='MYCa',loci=2,status='finishing',published=True)]))
            page=(q.root/'delivery/index.html').read_text()
            self.assertIn('MYCa',page)

    @unittest.skipUnless(Path('/proc').exists(),'Linux memory guard integration')
    def test_ramp_spaces_resumes_and_memory_guard_pauses_real_processes(self):
        children=[subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)']) for _ in range(2)]
        try:
            with TemporaryDirectory() as t:
                q=Queue(t)
                policy=dict(DEFAULTS,max_active=8,max_threads=256,active_growth_reserve_bytes=16*1024**3,pause_memory_bytes=128*1024**3,admission_interval_seconds=30)
                jobs=[]
                for i,child in enumerate(children):
                    ref=process(child.pid);pause_tree(ref)
                    jobs.append(dict(id='Q'+str(i),name='synthetic',panel='test',loci=0,status='paused',process=ref,threads=32,output=t+'/absent'))
                atomic(q.root/'state.json',dict(policy=policy,jobs=jobs,cgroup='synthetic',oom_kill_baseline=0))
                mem=dict(hard_limit=str(policy['memory_max_bytes']),soft_limit=str(policy['memory_high_bytes']),swap_limit='0',nonreclaimable_bytes=20*1024**3,events=dict(oom_kill=0))
                with patch('genecontext.stage2_queue.memory_state',return_value=mem):
                    q.step();q.step()
                    saved=json.loads((q.root/'state.json').read_text())
                    self.assertEqual([j['status'] for j in saved['jobs']],['running','paused'])
                    saved['last_admission']=time.time()-31;atomic(q.root/'state.json',saved)
                    q.step();saved=json.loads((q.root/'state.json').read_text())
                    self.assertEqual([j['status'] for j in saved['jobs']],['running','running'])
                    mem['nonreclaimable_bytes']=129*1024**3
                    q.step();saved=json.loads((q.root/'state.json').read_text())
                    self.assertTrue(saved['memory_hold'])
                    self.assertTrue(all(j['status']=='paused' for j in saved['jobs']))
                    self.assertTrue(all(process(c.pid)['state']=='T' for c in children))
        finally:
            for child in children:child.kill();child.wait()

    @unittest.skipUnless(Path('/proc').exists(),'Linux process adoption test')
    def test_pause_resume_preserves_process_identity_and_work_without_relaunch(self):
        child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(20)'])
        try:
            ref=process(child.pid);self.assertTrue(matches(ref))
            pause_tree(ref)
            deadline=time.time()+2
            while process(child.pid)['state']!='T' and time.time()<deadline:time.sleep(.01)
            self.assertEqual(process(child.pid)['state'],'T')
            resume_tree(ref);self.assertEqual(matches(ref)['start_ticks'],ref['start_ticks'])
            stale=dict(ref,start_ticks=ref['start_ticks']+1);self.assertIsNone(matches(stale))
            with self.assertRaisesRegex(RuntimeError,'identity changed'):resume_tree(stale)
        finally:child.kill();child.wait()

if __name__=='__main__':unittest.main()
