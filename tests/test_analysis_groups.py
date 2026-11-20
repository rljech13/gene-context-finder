"""Separate reviewed paralog sets without repeating discovery or mixing windows."""
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from genecontext.analysis_groups import normalize_groups
from genecontext.cli import parser
from genecontext.sequence import sha256, write_tsv
from genecontext.stage2 import selected_regions, plan


def fixture():
    loci=[dict(genome_id=g,locus_id=l,gene='FAMILY',contig='chr',strand='+',
               start=n,end=n+50,contig_length=1000,pass_filters=True,boundary_source='annotation')
          for g,l,n in [('g1','a1',100),('g1','b1',300),('g2','a2',100),('g2','c2',180),('g3','b3',300)]]
    data=dict(schema='blast_loci_v1',review_id='review',config=dict(gene='FAMILY',window=100),
              genomes=[dict(genome_id=g,status='completed') for g in ['g1','g2','g3']],loci=loci)
    groups=[dict(group_id=gid,name=name,members=[dict(genome_id=g,locus_id=l) for g,l in members])
            for gid,name,members in [('A','Family A',[('g1','a1'),('g2','a2')]),
                                     ('B','Family B',[('g1','b1'),('g3','b3')]),
                                     ('C','Family C',[('g2','c2')])]]
    draft=dict(schema='stage1_pair_draft_v2',review_id='review',analysis_groups=groups,
               pairs=[dict(genome_id=g,action='inspect',committed=g!='g3',
                           selected=[l['locus_id'] for l in loci if l['genome_id']==g],
                           reason='compare_copies',notes='Manual synteny grouping',overrides={}) for g in ['g1','g2','g3']],
               neighbor_shortlist=[dict(genome_id=g,locus_id=l) for g,l in [('g1','a1'),('g1','b1')]])
    return data,draft


class AnalysisGroups(unittest.TestCase):
    def args(self,path,**changes):
        return SimpleNamespace(selection=None,decisions=None,draft=str(path),committed_only=False,analysis_group=None,**changes)

    def test_requires_explicit_group_and_scoped_export_selects_it(self):
        data,draft=fixture();before=deepcopy((data,draft))
        with TemporaryDirectory() as temp:
            path=Path(temp)/'review.json';path.write_text(json.dumps(draft))
            args=self.args(path)
            with self.assertRaisesRegex(ValueError,'choose --analysis-group'):selected_regions(data,args)
            draft['active_analysis_group']='A';path.write_text(json.dumps(draft))
            regions,receipt,rows=selected_regions(data,args)
            self.assertEqual({r['locus_id'] for r in regions},{'a1','a2'})
            self.assertEqual({r['analysis_group'] for r in regions},{'Family A'})
            self.assertEqual(receipt['omitted_other_group_genomes'],['g3'])
            self.assertEqual(receipt['neighbor_shortlist'],[dict(genome_id='g1',locus_id='a1')])
            outside=[r for r in rows if r['locus_id'] in {'b1','c2'}]
            self.assertTrue(all(r['selection_status']=='outside_analysis_group' and r['original_decision']=='inspect' for r in outside))
            args.analysis_group='B'
            with self.assertRaisesRegex(ValueError,'differs'):selected_regions(data,args)
        self.assertEqual(data,before[0])

    def test_unsaved_members_never_approved_but_other_groups_may_be_pending(self):
        data,draft=fixture()
        with TemporaryDirectory() as temp:
            path=Path(temp)/'review.json';path.write_text(json.dumps(draft));args=self.args(path)
            args.analysis_group='Family A';self.assertEqual(len(selected_regions(data,args)[0]),2)
            args.analysis_group='B'
            with self.assertRaisesRegex(ValueError,'unsaved'):selected_regions(data,args)
            args.committed_only=True;regions,receipt,_=selected_regions(data,args)
            self.assertEqual([r['locus_id'] for r in regions],['b1'])
            self.assertEqual(receipt['omitted_unsaved_genomes'],['g3'])

    def test_unknown_unselected_or_multiply_assigned_members_rejected(self):
        data,draft=fixture();selected={(m['genome_id'],m['locus_id']) for g in draft['analysis_groups'] for m in g['members']}
        for mutation in ['unknown','wrong_genome','unselected','duplicate','path_id','duplicate_name']:
            groups=deepcopy(draft['analysis_groups']);eligible=set(selected)
            if mutation=='unknown':groups[0]['members'][0]['locus_id']='missing'
            elif mutation=='wrong_genome':groups[0]['members'][0]['genome_id']='g2'
            elif mutation=='unselected':eligible.remove(('g1','a1'))
            elif mutation=='duplicate':groups[1]['members'].append(groups[0]['members'][0])
            elif mutation=='path_id':groups[0]['group_id']='../escape'
            else:groups[1]['name']='family a'
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):normalize_groups(data,groups,eligible)

    def test_grouped_tsv_keeps_groups_and_checks_foreign_rows_outside_scope(self):
        data,draft=fixture();groups={m['locus_id']:g for g in draft['analysis_groups'] for m in g['members']}
        rows=[dict(review_id='review',genome_id=l['genome_id'],locus_id=l['locus_id'],decision='inspect',
                   analysis_group_id=groups[l['locus_id']]['group_id'],analysis_group=groups[l['locus_id']]['name']) for l in data['loci']]
        with TemporaryDirectory() as temp:
            path=Path(temp)/'review.tsv';write_tsv(path,rows,list(rows[0]))
            args=self.args(None);args.draft=None;args.decisions=str(path);args.analysis_group='A'
            self.assertEqual({r['locus_id'] for r in selected_regions(data,args)[0]},{'a1','a2'})
            rows[-1]['review_id']='another_checkpoint';write_tsv(path,rows,list(rows[0]))
            with self.assertRaisesRegex(ValueError,'checkpoint'):selected_regions(data,args)

    def test_plan_runs_each_group_separately_even_when_genomic_windows_overlap(self):
        data,draft=fixture();draft['pairs'][-1]['committed']=True
        with TemporaryDirectory() as temp:
            root=Path(temp);source=root/'stage1';source.mkdir()
            (source/'evidence.json').write_text(json.dumps(data))
            (source/'seal.json').write_text(json.dumps({'evidence.json':sha256(source/'evidence.json')}))
            path=root/'groups.json';path.write_text(json.dumps(draft))
            for gid,expected in [('A',{'a1','a2'}),('B',{'b1','b3'}),('C',{'c2'})]:
                args=parser().parse_args(['stage2-plan','--window-policy','pre-synteny','--results',str(source),'--draft',str(path),
                                         '--analysis-group',gid,'--out',str(root/gid)])
                self.assertEqual(plan(args),0)
                result=json.loads((root/gid/'window_plan.json').read_text())
                self.assertEqual({r['locus_id'] for r in result['regions']},expected)
                self.assertEqual(result['selection']['analysis_group']['group_id'],gid)
                self.assertTrue(all(set(p['region_ids'])<=expected for p in result['prediction_intervals']))

    def test_old_ungrouped_drafts_still_work(self):
        data,draft=fixture();draft['schema']='stage1_pair_draft_v1';draft.pop('analysis_groups')
        with TemporaryDirectory() as temp:
            path=Path(temp)/'old.json';path.write_text(json.dumps(draft));args=self.args(path);args.committed_only=True
            regions,receipt,_=selected_regions(data,args)
            self.assertEqual({r['locus_id'] for r in regions},{'a1','b1','a2','c2'})
            self.assertNotIn('analysis_group',receipt)
