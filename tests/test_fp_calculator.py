import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from src.fingerprint_calculator import FingerprintCalculator

class TestFingerprintCalculator(unittest.TestCase):
    def setUp(self):
        self.fp_calculator = FingerprintCalculator()
        self.test_smiles = [
            'CCCS(=O)(=O)N1CCC(CNC(=O)c2c(F)ccc(Cl)c2F)(C(=O)N2CCOCC2)CC1',
            'CCNC(=O)Nc1cn2c(-c3ncc(C)cn3)cc(-c3cncnc3)cc2n1',
            'COc1cccc(-c2cccc(-n3cc(CNCC4CCCCC4)c4ccccc43)c2)c1',
            'O=C(Cc1cc(Br)cc(Br)c1OC(=O)CCN1C(=O)CCC1=O)Nc1ccccc1[N+](=O)[O-]'
        ]

    def test_calculate_fingerprints(self):
        fingerprints = self.fp_calculator.calculate_fingerprints(self.test_smiles)
        
        # Check if the number of fingerprints matches the number of input SMILES
        self.assertEqual(len(fingerprints), len(self.test_smiles))
        
        # Check if identical SMILES produce identical fingerprints
        self.assertEqual(fingerprints[0], fingerprints[1])
        
        # Check if different SMILES produce different fingerprints
        self.assertNotEqual(fingerprints[0], fingerprints[2])
        self.assertNotEqual(fingerprints[0], fingerprints[3])
        self.assertNotEqual(fingerprints[2], fingerprints[3])

    def test_invalid_smiles(self):
        invalid_smiles = ['invalid_smiles_string', 'C1C']
        with self.assertRaises(ValueError):
            self.fp_calculator.calculate_fingerprints(invalid_smiles)

    def test_empty_input(self):
        empty_input = []
        fingerprints = self.fp_calculator.calculate_fingerprints(empty_input)
        self.assertEqual(len(fingerprints), 0)

    def test_fingerprint_consistency(self):
        # Test if the same SMILES always produces the same fingerprint
        smiles = 'CCNC(=O)Nc1cn2c(-c3ncc(C)cn3)cc(-c3cncnc3)cc2n1'
        fp1 = self.fp_calculator.calculate_fingerprints([smiles])[0]
        fp2 = self.fp_calculator.calculate_fingerprints([smiles])[0]
        self.assertEqual(fp1, fp2)

if __name__ == '__main__':
    unittest.main()