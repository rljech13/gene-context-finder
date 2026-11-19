"""3D presentation must not swap structures or inflate duplicate proteins."""
import hashlib,json,tempfile
from pathlib import Path
import unittest
from genecontext.stage4_viewer import write_viewer
from genecontext.stage4_structure import import_structures
from tests.test_stage4 import pdb,digest


class ViewerContract(unittest.TestCase):
    def test_exact_sequence_aliases_share_a_view_and_unknown_bfactor_is_not_plddt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sequence='M'+'A'*25;path=pdb(root,sequence)
            entry=dict(source='Synthetic test',file=path.name,file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),protein_sha256=digest(sequence))
            catalogue=root/'structures.json';catalogue.write_text(json.dumps(dict(schema='stage4_structures_v1',entries=[dict(entry,sequence_id=s) for s in ['model','reference_q']])))
            structures,_=import_structures(catalogue,{'model':sequence,'reference_q':sequence},root)
            data=dict(structures=structures,rows=[dict(model_id='model',species='Species',original_model='transcript')],trees=[dict(leaves=[dict(sequence_ids=['reference_q'],labels=['Supplied reference'])])])
            summary=write_viewer(root,data)
            self.assertEqual(summary['unique_structures'],1)
            self.assertEqual(summary['selected_models'],1)
            self.assertTrue((root/'structure_viewer.html').is_file())
            self.assertEqual(data['structures']['model']['confidence_kind'],'not_available')
            html=(root/'structure_viewer.html').read_text()
            self.assertIn('CA_trace',html);self.assertIn('Selected model',html)
            self.assertNotIn('src="https://',html)
            self.assertIn('3Dmol-LICENSE.txt',html)

    def test_changed_coordinates_rejected_before_viewer_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sequence='M'+'A'*25;path=pdb(root,sequence)
            entry=dict(source='Synthetic test',file=path.name,file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),protein_sha256=digest(sequence),sequence_id='model')
            cat=root/'structures.json';cat.write_text(json.dumps(dict(schema='stage4_structures_v1',entries=[entry])))
            structures,_=import_structures(cat,{'model':sequence},root)
            (root/structures['model']['full_file']).write_text('changed')
            data=dict(structures=structures,rows=[dict(model_id='model',species='Species',original_model='transcript')],trees=[])
            with self.assertRaisesRegex(ValueError,'coordinates changed'):write_viewer(root,data)


class StructuralColorContract(unittest.TestCase):
    def test_real_transform_and_residue_pairing_are_required(self):
        from genecontext.stage4_similarity import matrix_from_text,move,residue_distances
        matrix=matrix_from_text('0 2 0 -1 0\n1 3 1 0 0\n2 4 0 0 1\n')
        self.assertEqual(move([5,6,7],matrix),[-4,8,11])
        from genecontext.stage4_structure import ca_atoms
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);atoms=ca_atoms(pdb(root,'MA'))
            # Identity mapping with independently specified tool distances.
            identity=matrix_from_text('0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n')
            text='Aligned length= 2, RMSD= 0.00\nTM-score= 1.0 (normalized by length of Structure_1)\nTM-score= 1.0 (normalized by length of Structure_2)\n(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)\nMA\n::\nMA\n#Aligned atom 1\tAligned atom 2\tDistance#\nCA ALA A 1\tCA ALA A 1\t0.000\nCA MET A 2\tCA MET A 2\t0.000\n########\n'
            pairs,_=residue_distances(text,identity,atoms,atoms)
            self.assertEqual(len(pairs),2)
            self.assertEqual(pairs[0]['distance_angstrom'],0)
            with self.assertRaisesRegex(ValueError,'distances do not reproduce'):residue_distances(text,matrix,atoms,atoms)
            with self.assertRaisesRegex(ValueError,'sequence differs'):residue_distances(text.replace('\nMA\n','\nMM\n',1),identity,atoms,atoms)
            with self.assertRaisesRegex(ValueError,'reflection'):matrix_from_text('0 0 -1 0 0\n1 0 0 1 0\n2 0 0 0 1\n')

if __name__=='__main__':unittest.main()
