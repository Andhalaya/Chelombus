from multiprocessing import Pool
from typing import Optional
from config import N_JOBS
from rdkit import Chem
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import rdMolDescriptors                                                           
import numpy as np
import tmap as tm

class FingerprintCalculator:
    def __init__(self, smiles_list: list, fingerprint_type: str, permutations: Optional[int] = 512, fp_size: Optional[int] = 1024 ):
        self.smiles_list = smiles_list
        self.fingerprint_type = fingerprint_type.lower()
        self.permutations = permutations
        self.fp_size = fp_size

    def _calculate_mhfp_fp(self, smiles: str) -> np.array:
        """Calculate MHFP for a single SMILES string"""
        encoder = MHFPEncoder(self.permutations) 
        return tm.VectorUint(encoder.encode(smiles))
         
    def _calculate_mqn_fp(self, smiles: str) -> np.array: 
        """Calculate fp for a single SMILES string"""
        try:
            fingerprint = rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smiles))
            return np.array(fingerprint)
        
        except Exception as e: 
            print(f"error processing SMILES '{smiles}: {e}")
            return None
            
    def _calculate_mapc_fp(self, smiles: str) -> np.array:
        """Calculate MAP chiral fingerprints for a single SMILES string"""
        pass 
        
    def _fp_method(self):
        """
        Returns the correct function to calculate the fingerprints based on user selection
        """
        if self.fingerprint_type == 'mhfp':
            return self._calculate_mhfp_fp 
        elif self.fingerprint_type == 'mqn':
            return self._calculate_mqn_fp
        elif self.fingerprint_type == 'mapc':
            return self._calculate_mapc_fp
                    
    def calculate_fingerprints(self) -> np.array:
        fp_function = self._fp_method()
        with Pool(processes=N_JOBS) as pool: 
            fingerprints = pool.map(fp_function, self.smiles_list)
        return np.array(fingerprints)

"""
Example of usage: 
## Initiate 
fp_calculator = FingerprintCalculator(your_dataframe['smiles'], 'mqn') --> This will calculate the MQN fingerprints. For example is used when running the clustering pipeline
fp_calculator = FingerprintCalculator(smiles_list, 'mhfp', permutations=1024) --> This will calculate the MHFP with 1024 permutations (default is 512)

## Calculate fingerprints
fingerpints = fp_calculator.calculate_fingerprints() -> returns np.array witht the fp vectors of same length as smiles_list
"""
