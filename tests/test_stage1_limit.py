"""Per-assembly pool bounds, BLAST ordering, and immutable retained evidence."""
from collections import Counter
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from genecontext.cli import parser
from genecontext.stage1_limit import apply_limit, limit_candidates, subset_fasta


def pool(assemblies=18, copies=30):
    loci=[]
    for g in range(assemblies):
        for n in range(copies):
            loci.append(dict(locus_id=f'g{g}_{n}',genome_id=f'g{g}',contig='chr',start=n*100,
                rank_within_genome=n+1,chains=[dict(query='q',chain_id=f'c{g}_{n}',
                bitscore=(100000 if g==0 else 1000)-n,coverage=1,identity=1,evalue=1e-20)]))
    return dict(config=dict(gene='X'),queries=[dict(id='q',sequence='MA',members=[dict(query='Q',species='S')])],loci=loci)


class CandidateLimit(unittest.TestCase):
    def test_per_assembly_cap_retains_best_twenty_without_global_cutoff(self):
        data=pool(); original=deepcopy(data); audit=limit_candidates(data)
        self.assertEqual(len(data['loci']),360)
        counts=Counter(l['genome_id'] for l in data['loci'])
        self.assertEqual(len(counts),18)
        self.assertEqual(set(counts.values()),{20})
        for genome,count in counts.items():
            self.assertEqual({l['rank_within_genome'] for l in data['loci'] if l['genome_id']==genome},set(range(1,21)))
        by_id={l['locus_id']:l for l in original['loci']}
        self.assertTrue(all(l==by_id[l['locus_id']] for l in data['loci']))
        self.assertEqual(len(audit['omitted_loci']),180)
        reverse=deepcopy(original);reverse['loci'].reverse();limit_candidates(reverse)
        self.assertEqual({l['locus_id'] for l in data['loci']},{l['locus_id'] for l in reverse['loci']})

    def test_small_pools_are_unchanged_and_bound_is_hard(self):
        data=pool(18,5); before=deepcopy(data['loci']);limit_candidates(data)
        self.assertEqual(data['loci'],before)
        for invalid in (0,21,True,2.5):
            with self.assertRaises(ValueError):limit_candidates(pool(),invalid)
        with self.assertRaises(SystemExit):
            parser().parse_args(['stage1-limit','--results','a','--out','b','--max-loci-per-assembly','21'])

    def test_worker_defers_and_parent_bounds_the_combined_pool(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp); data=pool()
            for locus in data['loci']:
                locus.pop('rank_within_genome')  # Fresh search, before serialization.
            apply_limit(SimpleNamespace(_defer_candidate_limit=True),root,data)
            self.assertEqual(len(data['loci']),540)
            apply_limit(SimpleNamespace(),root,data)
            self.assertEqual(len(data['loci']),360)
            self.assertTrue((root/'raw/candidate_limit.json').is_file())

    def test_dna_export_contains_only_retained_records(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'loci.fna';path.write_text('>a detail\nAAA\n>b\nCCC\n')
            subset_fasta(path,{'b'})
            self.assertEqual(path.read_text(),'>b\nCCC\n')

    def test_limit_applies_before_early_synteny_when_disabled(self):
        from genecontext.stage1_context import attach
        with TemporaryDirectory() as tmp:
            data=pool();attach(SimpleNamespace(early_synteny='off'),Path(tmp),data)
            self.assertEqual(len(data['loci']),360)
