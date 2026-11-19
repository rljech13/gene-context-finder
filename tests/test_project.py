"""The project layer must preserve native evidence and explicit review gates."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from genecontext.cli import main, parser
from genecontext.project import Project, SCHEMA, attach, checkpoint, record_decision
from genecontext.project_run import plan, run_stage, status_rows
from genecontext.sequence import write_fasta, write_tsv, sha256
from test_stage3 import fixture


def seal(root):
    (root/'seal.json').write_text(json.dumps({str(p.relative_to(root)):sha256(p)
        for p in root.rglob('*') if p.is_file() and p.name!='seal.json'}))


def hashes(root):
    return {str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file()}


def study(root, gene='TEST'):
    work=root/'study';work.mkdir()
    write_fasta(root/'query.faa',[('q','MACDEFGHIKLMNPQRSTVWY')])
    write_fasta(root/'genome.fna',[('chr','ATGGCTTAA')])
    write_tsv(root/'genomes.tsv',[dict(genome_id='Anything_species',fasta='genome.fna')],['genome_id','fasta'])
    config=dict(schema=SCHEMA,name='Any study',genomes='../genomes.tsv',
        targets=[dict(id='custom',gene=gene,query='../query.faa',query_type='protein')],
        options=dict(stage1=dict(early_synteny='off',query_metadata_lookup='off'),stage4=dict(tree='none')))
    (work/'project.json').write_text(json.dumps(config))
    return Project(work)


def stage2(root):
    source=root/'native-stop2';source.mkdir()
    data,decisions,windows,coding,proteins=fixture()
    for name,records in [('windows.fna',windows),('cds.fna',coding),('proteins.faa',proteins)]:
        write_fasta(source/name,records.items())
    (source/'evidence.json').write_text(json.dumps(data))
    (source/'manifest.json').write_text(json.dumps(dict(stage='stage2',status='awaiting_manual_review')))
    seal(source)
    export=root/'decisions.json';export.write_text(json.dumps(decisions,ensure_ascii=False))
    return source,export


class ProjectSafety(unittest.TestCase):
    def test_config_is_generic_relative_and_status_plan_are_read_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root,'UNFAMILIAR_FAMILY');before=hashes(root)
            job=plan(p,'custom','stage1')
            self.assertIn('UNFAMILIAR_FAMILY',job['argv'])
            self.assertIn(str((root/'query.faa').resolve()),job['argv'])
            self.assertFalse(job['automatic_next_stage']);self.assertEqual(hashes(root),before)
            self.assertEqual(status_rows(p)[0]['status'],'ready')
            self.assertFalse(p.store.exists())
            with self.assertRaisesRegex(ValueError,'Complete or attach'):plan(p,'custom','stage3')

    def test_rejects_path_traversal_duplicate_ids_and_owned_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=study(Path(tmp));base=deepcopy(p.data)
            for change in ['traversal','duplicate','override','typo']:
                data=deepcopy(base)
                if change=='traversal':data['targets'][0]['id']='../escape'
                elif change=='duplicate':data['targets'].append(deepcopy(data['targets'][0]))
                elif change=='override':data['options']['stage1']['out']='../existing'
                else:data['options']['stage1']['threds']=8
                p.path.write_text(json.dumps(data))
                with self.assertRaises(ValueError):
                    other=Project(p.path);plan(other,'custom','stage1')

    def test_attach_and_decisions_preserve_all_original_bytes_and_need_matching_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root);source,export=stage2(root)
            originals=hashes(source);raw=export.read_bytes()
            run=attach(p,'custom','stage2',source)
            self.assertEqual(attach(p,'custom','stage2',source)['id'],run['id'])
            before=p.load();self.assertEqual(status_rows(p)[1]['status'],'awaiting_review')
            decision=record_decision(p,'custom','stage2',export)
            self.assertEqual((p.root/decision['object']).read_bytes(),raw)
            self.assertEqual(status_rows(p)[1]['status'],'reviewed')
            self.assertEqual(hashes(source),originals);self.assertEqual(export.read_bytes(),raw)
            bad=json.loads(raw);bad['review_id']='another-review';export.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError,'another'):record_decision(p,'custom','stage2',export)
            self.assertEqual(len(p.load()['decisions']),1)
            self.assertEqual(len(before['runs']),1)
            (source/'proteins.faa').write_text('>changed\nMA\n')
            self.assertEqual(status_rows(p)[1]['status'],'reviewed')
            # Full scientific verification at execution rejects changed evidence
            # even when the seal document itself was not edited.
            with self.assertRaisesRegex(ValueError,'changed'):plan(p,'custom','stage3')

    def test_no_silent_replacement_and_no_write_inside_imported_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root);source,_=stage2(root);attach(p,'custom','stage2',source)
            other=root/'other';shutil.copytree(source,other)
            data=json.loads((other/'evidence.json').read_text());data['review_id']='different'
            (other/'evidence.json').write_text(json.dumps(data));seal(other)
            with self.assertRaisesRegex(ValueError,'already exists'):attach(p,'custom','stage2',other)
            attach(p,'custom','stage2',other,True)
            self.assertEqual(len(p.load()['runs']),2)
            config=p.data;(source/'project.json').write_text(json.dumps(config))
            nested=Project(source)
            with self.assertRaisesRegex(ValueError,'outside'):attach(nested,'custom','stage2',source)

    def test_run_failure_keeps_attempt_and_never_retries_implicitly(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=study(Path(tmp))
            with patch('genecontext.project_run.subprocess.run') as execute:
                execute.return_value.returncode=9
                result=run_stage(p,'custom','stage1')
                self.assertEqual(result['status'],'failed');self.assertEqual(execute.call_count,1)
                with self.assertRaisesRegex(ValueError,'already exists'):run_stage(p,'custom','stage1')
                self.assertEqual(execute.call_count,1)
            self.assertEqual(len(p.load()['runs']),1)

    def test_html_escapes_names_and_does_not_start_calculations(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=study(Path(tmp),'<script>anything</script>')
            from genecontext.project_report import write_report
            with patch('subprocess.run',side_effect=AssertionError('No tools')):
                result=write_report(p)
            html=Path(result['report']).read_text();self.assertIn('&lt;script&gt;',html)
            self.assertFalse(p.store.exists())

    def test_changed_dataset_is_detected_without_replacing_saved_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root)
            first=plan(p,'custom','stage1')
            (root/'genome.fna').write_text('>chr\nATGGCTTAG\n')
            second=plan(p,'custom','stage1')
            self.assertNotEqual(first['config_sha256'],second['config_sha256'])

    def test_portable_package_contains_only_requested_platform_and_preserves_source(self):
        from genecontext.interactive_review import stage2_dataset,write_dataset,package_review
        import zipfile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,_=stage2(root)
            write_dataset(source,stage2_dataset(source,json.loads((source/'evidence.json').read_text())))
            seal(source);before=hashes(source)
            bins=root/'bins';bins.mkdir()
            for name in ['review-windows-arm64.exe','review-linux-amd64']:(bins/name).write_bytes(b'fixture launcher; never executed')
            (bins/'LICENSE-Go.txt').write_text('Fixture license')
            args=parser().parse_args(['review-package','--results',str(source),'--out',str(root/'app'),'--platform','windows-arm64'])
            with patch.dict('os.environ',GCF_REVIEW_BIN_DIR=str(bins)):
                package_review(args)
            with zipfile.ZipFile(root/'app/review_app.zip') as archive:
                binaries=[n for n in archive.namelist() if '/.review/bin/' in n]
                self.assertEqual(len(binaries),1);self.assertTrue(binaries[0].endswith('review-windows-arm64.exe'))
                self.assertIsNone(archive.testzip())
                saved=json.loads(archive.read('Gene_context_review/review_data.json'))
                self.assertTrue(saved['portable_available'])
            self.assertEqual(hashes(source),before)

    def test_shared_queue_reuses_existing_ticket_without_starting_a_dispatcher(self):
        from genecontext.project import atomic
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root);source,decision=stage2(root)
            queue=root/'queue';queue.mkdir()
            atomic(queue/'state.json',dict(policy=dict(max_job_threads=32),jobs=[]))
            repo=Path(__file__).resolve().parents[1]
            p.data['execution']=dict(stage2=dict(mode='queue',root=str(queue),source=str(repo),
                entrypoint=str(repo/'scripts/enqueue_stage2.py')))
            p.path.write_text(json.dumps(p.data));p=Project(p.path)
            prepared=dict(action='run',target='custom',stage='stage2',output=str(root/'future-output'),
                argv=[__import__('sys').executable,'-m','genecontext','stage2','--results',str(source),
                      '--decisions',str(decision),'--genomes',str(root/'genomes.tsv'),
                      '--out',str(root/'future-output'),'--threads','2'],
                config_sha256='a'*64,source_run_id=None,decision_sha256=sha256(decision),
                windows=dict(source_review_id='parent',regions=[dict(locus_id='locus')]),executor='queue')
            with patch('genecontext.project_run.plan',return_value=deepcopy(prepared)):
                first=run_stage(p,'custom','stage2')
                prepared['output']=str(root/'deduplicated-output')
                prepared['argv'][prepared['argv'].index('--out')+1]=prepared['output']
            with patch('genecontext.project_run.plan',return_value=deepcopy(prepared)):
                second=run_stage(p,'custom','stage2',True)
            self.assertEqual(first['queue_receipt']['id'],second['queue_receipt']['id'])
            self.assertTrue(first['queue_receipt']['created']);self.assertFalse(second['queue_receipt']['created'])
            self.assertEqual(first['output'],second['output'])
            state=json.loads((queue/'state.json').read_text())
            self.assertEqual(len(state['jobs']),1);self.assertEqual(state['jobs'][0]['status'],'queued')
            self.assertFalse((root/'future-output').exists())


@unittest.skipUnless(shutil.which('mafft'),'requires actual MAFFT')
class ProjectHandoff(unittest.TestCase):
    def test_real_comparison_and_final_dossier_preserve_inputs_and_review_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root);source,export=stage2(root);original=hashes(source)
            attach(p,'custom','stage2',source)
            with self.assertRaisesRegex(ValueError,'explicit'):plan(p,'custom','stage3')
            record_decision(p,'custom','stage2',export)
            run3=run_stage(p,'custom','stage3');self.assertEqual(run3['status'],'awaiting_review',Path(run3['log']).read_text())
            self.assertFalse((Path(run3['output'])/'review_app.zip').exists())
            self.assertEqual(plan(p,'custom','stage3')['action'],'reuse')
            with self.assertRaisesRegex(ValueError,'explicit'):plan(p,'custom','stage4')
            final=json.loads((Path(run3['output'])/'finalization.decisions.json').read_text())
            for choice in final['choices'].values():choice['reviewed']=True;choice['reason']='Explicit software control review'
            choice=root/'final.json';choice.write_text(json.dumps(final))
            record_decision(p,'custom','stage3',choice)
            run4=run_stage(p,'custom','stage4');self.assertEqual(run4['status'],'completed',Path(run4['log']).read_text())
            self.assertEqual(hashes(source),original)
            # A new review version invalidates downstream display state but
            # never edits or removes any former output or decision object.
            old3=hashes(Path(run3['output']));old4=hashes(Path(run4['output']))
            changed=json.loads(export.read_text());next(iter(changed['regions'].values()))['notes']='A revised explicit comment'
            export.write_text(json.dumps(changed));record_decision(p,'custom','stage2',export)
            rows=status_rows(p);self.assertEqual(rows[2]['status'],'stale');self.assertEqual(rows[3]['status'],'stale')
            self.assertEqual(hashes(Path(run3['output'])),old3);self.assertEqual(hashes(Path(run4['output'])),old4)
            with self.assertRaisesRegex(ValueError,'already exists'):plan(p,'custom','stage3')
            new=run_stage(p,'custom','stage3',True)
            self.assertNotEqual(new['output'],run3['output']);self.assertEqual(len(p.load()['runs']),4)


@unittest.skipUnless(all(shutil.which(v) for v in ['samtools','makeblastdb','blastp']),'requires actual BLAST+ and samtools')
class ArbitraryInputs(unittest.TestCase):
    def test_two_unrelated_genes_two_species_and_stop1_gate(self):
        import random
        from genecontext.sequence import CODE
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=study(root);rng=random.Random(20261119)
            aa='ACDEFGHIKLMNPQRSTVWY';proteins=['M'+''.join(rng.choice(aa) for _ in range(99)) for _ in range(2)]
            codons={a:next(c for c,value in CODE.items() if value==a) for a in aa}
            positions=[];dna='A'*100
            for i,protein in enumerate(proteins):
                start=len(dna);dna+=''.join(codons[a] for a in protein)+'TAA';positions.append((start,len(dna)));dna+='A'*100
            rows=[]
            for species in ['Species_one','Species_two']:
                write_fasta(root/(species+'.fna'),[('chr',dna)])
                lines=['##gff-version 3',f'##sequence-region chr 1 {len(dna)}']
                for i,(start,end) in enumerate(positions):
                    lines += [f'chr\ttest\tgene\t{start+1}\t{end}\t.\t+\t.\tID=g{i};Name=GENE{i}',f'chr\ttest\tmRNA\t{start+1}\t{end}\t.\t+\t.\tID=t{i};Parent=g{i}',f'chr\ttest\tCDS\t{start+1}\t{end}\t.\t+\t0\tParent=t{i}']
                (root/(species+'.gff3')).write_text('\n'.join(lines)+'\n')
                write_fasta(root/(species+'.faa'),[(f't{i}',v) for i,v in enumerate(proteins)])
                rows.append(dict(genome_id=species,fasta=species+'.fna',gff3=species+'.gff3',protein=species+'.faa'))
            write_tsv(root/'genomes.tsv',rows,list(rows[0]))
            p.data['targets']=[]
            for i,name in enumerate(['GENE_ALPHA','GENE_BETA']):
                write_fasta(root/(name+'.faa'),[('q',proteins[i])])
                p.data['targets'].append(dict(id=name.lower(),gene=name,query='../'+name+'.faa',query_type='protein'))
            p.path.write_text(json.dumps(p.data));p=Project(p.path)
            before=hashes(root);before={k:v for k,v in before.items() if not k.startswith('study/.gcf')}
            for tid in p.targets:
                result=run_stage(p,tid,'stage1');self.assertEqual(result['status'],'awaiting_review',Path(result['log']).read_text())
                evidence=json.loads((Path(result['output'])/'evidence.json').read_text())
                self.assertEqual({l['genome_id'] for l in evidence['loci']},{'Species_one','Species_two'})
                with self.assertRaisesRegex(ValueError,'explicit'):plan(p,tid,'stage2')
            after={k:v for k,v in hashes(root).items() if not k.startswith('study/.gcf')}
            self.assertEqual(before,after)


if __name__=='__main__':unittest.main()
