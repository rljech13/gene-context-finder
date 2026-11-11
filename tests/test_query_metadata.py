"""Source context must not depend on manually reconstructed per-gene metadata."""
from pathlib import Path
from copy import deepcopy
import hashlib
import tempfile
import unittest
from genecontext.query_metadata import fasta_metadata, restore_headers
from genecontext.cli import parser
from genecontext.stage1 import panel
from genecontext.early_synteny import resolve_origins
from genecontext.sequence import sha256


class QuerySources(unittest.TestCase):
    def test_standard_headers_are_retained_and_explicit_tsv_takes_precedence(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);f=root/'q.faa'
            f.write_text('>sample_NP_123.1 protein [organism=Species one] [GeneID=123]\nM'+'A'*40+'\n>other protein [Species two]\nM'+'C'*40+'\n')
            info=fasta_metadata(f);self.assertEqual(info['sample_NP_123.1']['species'],'Species one');self.assertEqual(info['sample_NP_123.1']['gene_id'],'123');self.assertEqual(info['sample_NP_123.1']['protein'],'NP_123.1')
            args=parser().parse_args(['stage1','--gene','GENE','--query',str(f),'--query-type','protein','--genomes','unused','--out','unused'])
            queries,_=panel(args);self.assertEqual(queries[0]['members'][0]['gene_id'],'123');self.assertIn('[organism=',queries[0]['members'][0]['source_header'])
            meta=root/'meta.tsv';meta.write_text('query\tspecies\tgene_id\nsample_NP_123.1\tSpecies curated\t456\nother\tSpecies two\t\n');args.query_metadata=str(meta)
            queries,_=panel(args);self.assertEqual(queries[0]['members'][0]['species'],'Species curated');self.assertEqual(queries[0]['members'][0]['gene_id'],'456')

    def test_header_migration_preserves_query_sequences_and_search_chains(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'q.faa';p.write_text('>a protein [Species one] [GeneID=1]\nM'+'A'*40+'\n')
            data=dict(query_provenance=dict(sha256=sha256(p)),queries=[dict(id='q1',sequence='M'+'A'*40,members=[dict(query='a')])],loci=[dict(query_ids=['q1'],chains=[dict(bitscore=50,hsps=[dict(qstart=0,qend=40)])])]);chains=deepcopy(data['loci'][0]['chains'])
            restore_headers(data,p);self.assertEqual(data['queries'][0]['members'][0]['gene_id'],'1');self.assertEqual(data['loci'][0]['support_species'],1);self.assertEqual(data['loci'][0]['chains'],chains)
            p.write_text(p.read_text().replace('GeneID=1','GeneID=2'))
            with self.assertRaisesRegex(ValueError,'differs'):restore_headers(data,p)

    def test_ambiguous_metadata_is_not_guessed(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'q.faa';p.write_text('>q [organism=Species one] [organism=Species two] [GeneID=1] [GeneID=2]\nMAAA\n')
            value=fasta_metadata(p)['q'];self.assertNotIn('species',value);self.assertNotIn('gene_id',value)

    def catalog(self,gid='a',gene='gene1'):
        return dict(genome_id=gid,assembly=gid,species='Species one',aliases={'NP_123.1':[(gene,'old_sequence_hash')]},sequence_genes={},ncbi_gene_ids={'123':[gene]})

    def test_unique_exact_accession_binds_context_but_never_verifies_different_sequence(self):
        q=dict(id='q',sequence='MAAAA',members=[dict(query='sample_NP_123.1',species='Species one')]);refs,rows=resolve_origins([q],[self.catalog()])
        self.assertEqual(len(refs),1);self.assertEqual(rows[0]['status'],'identifier_linked_sequence_differs');self.assertTrue(rows[0]['context_binding_only']);self.assertFalse(rows[0]['query_sequence_verified']);self.assertEqual(rows[0]['binding_basis'],'exact_versioned_accession')
        refs,rows=resolve_origins([q],[self.catalog(),self.catalog('b','gene2')]);self.assertEqual(refs,[]);self.assertFalse(rows[0]['context_binding_only']);self.assertEqual(len(rows[0]['mismatching_sources']),2)
        q['members'][0]['species']='Species other';self.assertEqual(resolve_origins([q],[self.catalog()])[0],[])

    def test_unversioned_symbol_or_unknown_identifier_cannot_bind_a_mismatch(self):
        catalog=self.catalog();catalog['aliases']['GENE']=catalog['aliases']['NP_123.1'];q=dict(id='q',sequence='MAAAA',members=[dict(query='GENE')]);self.assertEqual(resolve_origins([q],[catalog])[0],[])

    def test_explicit_header_gene_id_binds_across_isoform_disagreement(self):
        q=dict(id='q',sequence='MAAAA',members=[dict(query='different',species='Species one',gene_id='123')]);refs,rows=resolve_origins([q],[self.catalog()]);self.assertEqual(refs[0]['gene_id'],'gene1');self.assertEqual(rows[0]['binding_basis'],'explicit_gene_id');self.assertFalse(rows[0]['query_sequence_verified'])


class RegisteredSources(unittest.TestCase):
    def test_shared_inventory_is_automatic_and_explicit_input_overrides_it(self):
        from types import SimpleNamespace
        from unittest.mock import patch
        from genecontext.stage1_context import attach
        for kind in ['registered','explicit','none']:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as t:
                root=Path(t);out=root/'out';(out/'raw').mkdir(parents=True);cache=root/'cache';cache.mkdir();targets=root/'targets.tsv';targets.write_text('genome_id\tfasta\n');explicit=root/'explicit.tsv';explicit.write_text('genome_id\tfasta\n')
                if kind!='none':(cache/'source_datasets.tsv').write_text('genome_id\tfasta\n')
                args=SimpleNamespace(genomes=str(targets),neighbor_cache=str(cache),source_datasets=str(explicit) if kind=='explicit' else None,threads=1)
                data=dict(review_id='review',config=dict(query_type='dna'),query_provenance={},queries=[],genomes=[],loci=[])
                result=dict(schema='early_synteny_v1',loci=[],errors=[],source_resolution=[])
                with patch('genecontext.early_synteny.run_early_synteny',return_value=result) as run:
                    self.assertFalse(attach(args,out,data))
                selected=run.call_args.kwargs['source_genomes_manifest'];expected=str(explicit) if kind=='explicit' else cache/'source_datasets.tsv' if kind=='registered' else None
                self.assertEqual(selected,expected)
