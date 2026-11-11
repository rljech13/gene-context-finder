import hashlib,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from io import BytesIO
from genecontext.query_source_lookup import enrich

XML=b'''<GBSet><GBSeq><GBSeq_accession-version>NP_123.1</GBSeq_accession-version><GBSeq_organism>Species one</GBSeq_organism><GBSeq_sequence>MAAAA</GBSeq_sequence><GBSeq_feature-table><GBFeature><GBFeature_quals><GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>GeneID:123</GBQualifier_value></GBQualifier><GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>GENE</GBQualifier_value></GBQualifier></GBFeature_quals></GBFeature></GBSeq_feature-table></GBSeq></GBSet>'''

class RefSeqSources(unittest.TestCase):
 def test_exact_sequence_metadata_is_cached_and_usable_offline(self):
  with tempfile.TemporaryDirectory() as d:
   r=Path(d);queries=[dict(sequence='MAAAA',members=[dict(query='sample_NP_123.1')])]
   with patch('genecontext.query_source_lookup.urlopen',return_value=BytesIO(XML)),patch('genecontext.query_source_lookup.time.sleep'):
    result=enrich(queries,r/'cache',r/'out')
   self.assertEqual(queries[0]['members'][0]['gene_id'],'123');self.assertEqual(queries[0]['members'][0]['species'],'Species one');self.assertEqual(queries[0]['sequence'],'MAAAA')
   fresh=[dict(sequence='MAAAA',members=[dict(query='sample_NP_123.1')])]
   with patch('genecontext.query_source_lookup.urlopen',side_effect=AssertionError('No network')):enrich(fresh,r/'cache',r/'offline','off')
   self.assertEqual(fresh[0]['members'][0]['gene_id'],'123')
   (r/'cache/NP_123.1.xml').write_bytes(b'corrupt')
   result=enrich([dict(sequence='MAAAA',members=[dict(query='sample_NP_123.1')])],r/'cache',r/'bad')
   self.assertEqual(result['records'][0]['status'],'metadata_lookup_unavailable')
 def test_wrong_sequence_species_and_accession_never_supply_gene_id(self):
  for changes,raw in [({'sequence':'MCCCC'},XML),({'species':'Species other'},XML),({},XML.replace(b'NP_123.1',b'NP_999.1'))]:
   with self.subTest(changes=changes),tempfile.TemporaryDirectory() as d:
    r=Path(d);m=dict(query='q_NP_123.1');m.update({k:v for k,v in changes.items() if k!='sequence'});q=dict(sequence=changes.get('sequence','MAAAA'),members=[m])
    with patch('genecontext.query_source_lookup.urlopen',return_value=BytesIO(raw)),patch('genecontext.query_source_lookup.time.sleep'):result=enrich([q],r/'cache',r/'out')
    self.assertNotIn('gene_id',m);self.assertNotEqual(result['records'][0]['status'],'verified_refseq_metadata')
 def test_off_mode_does_not_contact_ncbi_or_invent_a_source(self):
  with tempfile.TemporaryDirectory() as d:
   r=Path(d);q=dict(sequence='MAAAA',members=[dict(query='q_NP_123.1')])
   with patch('genecontext.query_source_lookup.urlopen',side_effect=AssertionError('No network')):result=enrich([q],r/'cache',r/'out','off')
   self.assertEqual(result['records'][0]['status'],'lookup_disabled_cache_missing');self.assertNotIn('gene_id',q['members'][0])
