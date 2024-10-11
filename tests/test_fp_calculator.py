import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import unittest
from unittest.mock import patch, MagicMock
from multiprocessing import Pool
from src.fingerprint_calculator import FingerprintCalculator, _calculate_fingerprint, init_worker
from config import N_JOBS

class TestFingerprintCalculator(unittest.TestCase):

    def setUp(self):
        self.calculator = FingerprintCalculator()
        self.sample_smiles = [
            'CCNC(=O)Nc1cn2c(-c3ncc(C)cn3)cc(-c3cncnc3)cc2n1', 
            'CN1C(=O)[C@@H]2CN(S(=O)(=O)c3ccc4occc4c3)C[C@@H]2c2cnccc21', 
            'COC(=O)[C@@H]1[C@@H]2C(=O)N(CC(=O)NC[C@H]3O[C@@H](n4ccc(=O)[nH]c4=O)[C@H](O)[C@@H]3O)C(=O)[C@@H]2N[C@@H]1c1ccccc1',
            'C#Cc1ccc(-c2ccccc2)cc1',
            'Cc1cccc(C(=O)NC2(C(=O)O)Cc3ccc(F)cc3C2)c1C1=CCCC1']

    @patch('src.fingerprint_calculator.Pool')
    def test_calculate_fingerprints(self, mock_pool):
        # Mock the Pool.map method to return some dummy fingerprints
        mock_pool.return_value.__enter__.return_value.map.return_value = [
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [1, 1, 0, 0],
            [0, 0, 1, 1]
        ]

        result = self.calculator.calculate_fingerprints(self.sample_smiles)

        # Check if the result is as expected
        self.assertEqual(len(result), len(self.sample_smiles))
        self.assertTrue(all(isinstance(fp, list) for fp in result))

        # Check if Pool was called with the correct arguments
        mock_pool.assert_called_once_with(processes=N_JOBS, initializer=init_worker)
        mock_pool.return_value.__enter__.return_value.map.assert_called_once_with(_calculate_fingerprint, self.sample_smiles)

    @patch('src.fingerprint_calculator.MHFPEncoder')
    def test_calculate_fingerprint_valid_smiles(self, mock_encoder_class):
        # Mock the MHFPEncoder instance
        mock_encoder = MagicMock()
        mock_encoder.encode.return_value = [1, 0, 1, 0]
        mock_encoder_class.return_value = mock_encoder

        # Call init_worker to set up the global encoder
        init_worker()

        result = _calculate_fingerprint('C')

        self.assertEqual(result, [1, 0, 1, 0])
        mock_encoder.encode.assert_called_once_with('C')

    @patch('src.fingerprint_calculator.MHFPEncoder')
    @patch('src.fingerprint_calculator.print')  # Mock the print function
    def test_calculate_fingerprint_invalid_smiles(self, mock_print, mock_encoder_class):
        # Mock the MHFPEncoder instance
        mock_encoder = MagicMock()
        mock_encoder.encode.side_effect = Exception("Invalid SMILES")
        mock_encoder_class.return_value = mock_encoder

        # Call init_worker to set up the global encoder
        init_worker()

        result = _calculate_fingerprint('Invalid SMILES')

        self.assertIsNone(result)
        mock_print.assert_called_once_with("error processing SMILES 'Invalid SMILES: Invalid SMILES")

if __name__ == '__main__':
    unittest.main()