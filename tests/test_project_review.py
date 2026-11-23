"""The current study route must be generic, reproducible and source-preserving."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from genecontext.project import Project,attach
from genecontext.project_review import prepare,record_handoff,handoff
from genecontext.project_run import plan,run_stage,status_rows
from genecontext.legacy_review import import_analysis
from genecontext.stage2_report import write_report
from test_project import study,stage2,seal,hashes


def workspace(root):
    project=study(root)
    data=project.data;data['workflow']='sequential_v2';data['options']['stage4']={};project.path.write_text(json.dumps(data));project=Project(project.path)
    source,choices=stage2(root);native=json.loads((source/'evidence.json').read_text())
    write_report(source,native);seal(source)
    attach(project,'custom','stage2',source)
    review=prepare(project,'custom');return project,source,choices,review


def exported(project,review,choices,path):
    pool=json.loads((project.root/review['directory']/'review_data.json').read_text());members=pool['sequences']
    decision=json.loads(choices.read_text());group=members[0]['group']
    protein={s['id']:s['protein'] for s in members};analyses=import_analysis(group,{'protein':protein},{s['id']:s for s in members},pool['review_id'],{})
    native=json.loads((Path(project.current(project.load(),'custom','stage2')['output'])/'evidence.json').read_text())
    data=dict(schema='gcf_stop4_input_v1',source_review_id=pool['source_review_id'],alignment_review_id=pool['review_id'],
        stage2=decision,stage3=dict(reviewed=True,decisions={s['id']:dict(role='primary_target',reason='Temporary software-control choice') for s in members},
        synteny=dict(regions=native['regions'],neighborhood=[dict(model_id=s['id'],region_id=s['metadata']['region_id'],physical_id=s['metadata']['model_constraints']['physical_id'],**{k:s['metadata'][k] for k in ('source','contig','start','end','strand')}) for s in members],reviewed_pairs={})),members=members,analyses={group:analyses},test_fixture_only=True)
    save_export(path,data);return data


def save_export(path,data):
    text=json.dumps(data);path.write_text(json.dumps(dict(payload=text,sha256=hashlib.sha256(text.encode()).hexdigest())))


class SequentialStudy(unittest.TestCase):
    def test_result_comments_are_preserved_and_malformed_comments_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,_,choices,review=workspace(root);path=root/'export.json'
            data=exported(p,review,choices,path)
            group=next(iter(data['analyses']))
            data['analyses'][group]['comments']={'alignment':'Inspect <gap> & CDS','tree':'Tree note'}
            save_export(path,data);record_handoff(p,'custom',path)
            result=run_stage(p,'custom','stage4')
            final=json.loads((Path(result['output'])/'evidence.json').read_text())
            self.assertEqual(final['trees'][0]['review_comments'],data['analyses'][group]['comments'])
            data['analyses'][group]['comments']={'tree':42}
            save_export(path,data)
            with self.assertRaisesRegex(ValueError,'review comment'):record_handoff(p,'custom',path)

    def test_prepared_models_to_final_dossier_without_recomputation_or_source_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,source,choices,review=workspace(root);before=hashes(source)
            self.assertEqual(len(p.load()['reviews']),1)
            self.assertEqual(prepare(p,'custom')['id'],review['id'])
            self.assertFalse((p.root/review['directory']/'review_app.zip').exists())
            page=(p.root/review['directory']/'report.html').read_text()
            self.assertIn('gcf-separate-synteny-stop-v2',page)
            self.assertIn('gcf-synteny-reference-approval-v1',page)
            self.assertNotIn('DNMT3B',page)
            with self.assertRaisesRegex(ValueError,'project handoff'):plan(p,'custom','stage4')
            export=root/'export.json';exported(p,review,choices,export)
            entry=record_handoff(p,'custom',export)
            self.assertEqual((p.root/entry['object']).read_bytes(),export.read_bytes())
            result=run_stage(p,'custom','stage4');self.assertEqual(result['status'],'completed',result)
            self.assertEqual(plan(p,'custom','stage4')['action'],'reuse')
            output=Path(result['output']);evidence=json.loads((output/'evidence.json').read_text())
            self.assertEqual(len(evidence['rows']),2)
            self.assertEqual(json.loads((output/'manifest.json').read_text())['commands'],[])
            self.assertEqual(hashes(source),before)
            self.assertEqual([r for r in status_rows(p) if r['stage']=='stage4'][0]['status'],'completed')

    def test_foreign_workspace_and_self_consistent_modified_sequence_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,_,choices,review=workspace(root);path=root/'export.json';base=exported(p,review,choices,path)
            data=deepcopy(base);data['alignment_review_id']='different';save_export(path,data)
            with self.assertRaises(ValueError):record_handoff(p,'custom',path)
            data=deepcopy(base);member=data['members'][0]
            member['cds']=member['cds'][:3]+'GTT'+member['cds'][6:]
            member['protein']=member['protein'][:1]+'V'+member['protein'][2:]
            for key in ('cds','protein'):member['metadata'][key+'_sha256']=hashlib.sha256(member[key].encode()).hexdigest()
            group=member['group'];data['analyses'][group]=import_analysis(group,{'protein':{s['id']:s['protein'] for s in data['members']}},
                {s['id']:s for s in data['members']},data['alignment_review_id'],{})
            save_export(path,data)
            from genecontext.stage4_handoff import read_handoff
            read_handoff(path)  # Internally consistent, but not the project's exact native models.
            with self.assertRaisesRegex(ValueError,'differs from the prepared workspace'):record_handoff(p,'custom',path)
            self.assertEqual(p.load().get('handoffs',[]),[])

    def test_source_context_tampering_and_changed_workspace_are_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,_,choices,review=workspace(root);path=root/'export.json';base=exported(p,review,choices,path)
            data=deepcopy(base);data['stage3']['synteny']['regions'][0]['anchor_start']+=1;save_export(path,data)
            with self.assertRaisesRegex(ValueError,'Synteny region differs'):record_handoff(p,'custom',path)
            save_export(path,base);record_handoff(p,'custom',path)
            old_sha=handoff(p,p.load(),'custom')['sha256']
            next_review=prepare(p,'custom',True);self.assertNotEqual(next_review['id'],review['id'])
            with self.assertRaisesRegex(ValueError,'earlier review'):plan(p,'custom','stage4')
            self.assertTrue(any(h['sha256']==old_sha for h in p.load()['handoffs']))

    def test_comment_escaping_and_synteny_overview_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            from genecontext.project_report import write_report as overview
            from genecontext.legacy_review import read_payload
            root=Path(tmp);p,source,_,review=workspace(root)
            raw=json.loads((source/'evidence.json').read_text());raw['regions'][0]['notes']='<gene> & </script> comment'
            # Re-render a new source fixture rather than editing an imported one.
            other=root/'new-source';import shutil;shutil.copytree(source,other)
            (other/'evidence.json').write_text(json.dumps(raw));write_report(other,raw);seal(other)
            attach(p,'custom','stage2',other,True);review=prepare(p,'custom',True)
            rendered=read_payload(p.root/review['directory']/'report.html')
            self.assertEqual(rendered['regions'][0]['notes'],raw['regions'][0]['notes'])
            summary=Path(overview(p)['report']).read_text()
            self.assertIn('?view=synteny',summary)

    def test_new_route_does_not_execute_the_legacy_stage3(self):
        with tempfile.TemporaryDirectory() as tmp:
            p,_,_,_=workspace(Path(tmp))
            with self.assertRaisesRegex(ValueError,'synteny review page'):plan(p,'custom','stage3')

if __name__=='__main__':unittest.main()
