from multiprocessing import Pool
from typing import Optional
from config import N_JOBS
from rdkit import Chem
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors                                                           
import numpy as np
import tmap as tm
import hashlib 
from mapchiral.mapchiral import encode as mapc_enc
import numpy.typing as npt

class FingerprintCalculator:
    def __init__(
            self,
            smiles_list: list, 
            fingerprint_type: str, 
            permutations: Optional[int] = 512, 
            fp_size: Optional[int] = 1024,
            radius: Optional[int] = 2):

        self.smiles_list = smiles_list
        self.fingerprint_type = fingerprint_type.lower()
        self.permutations = permutations
        self.fp_size = fp_size
        self.radius = radius

    def _calculate_mhfp_fp(self, smiles: str) -> np.array:
        """Calculate MHFP for a single SMILES string"""
        encoder = MHFPEncoder(self.permutations) 
        #TODO: Fix error. Can't picke tmap.VectorUint objects.
        return np.array(encoder.encode(smiles))
         
    def _calculate_mqn_fp(self, smiles: str) -> np.array: 
        """Calculate fp for a single SMILES string"""
        try:
            fingerprint = rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smiles))
            return np.array(fingerprint)
        
        except Exception as e: 
            print(f"error processing SMILES '{smiles}: {e}")
            return None   
            
    def _calculate_morgan_fp(self, smiles: str) -> np.array:
        """
        Calculate Morgan-Hashed Fingerprint (MHFP) for a single SMILES string.

        Args:
            smiles (str): SMILES representation of the molecule.
            radius (int): Radius for Morgan fingerprint.
            nBits (int): Number of bits in the fingerprint.

        Returns:
            np.array: MHFP as a numpy array.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string.")

            # Generate standard Morgan fingerprint
            morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=self.radius, nBits=self.fp_size)

            # Convert to bit string
            bit_string = morgan_fp.ToBitString()

            # Hash the bit string using SHA256 and truncate to nBits
            hash_obj = hashlib.sha256(bit_string.encode())
            hashed_bits = bin(int(hash_obj.hexdigest(), 16))[2:].zfill(256)[:self.fp_size]

            # Convert hashed bits to numpy array
            fingerprint = np.array([int(bit) for bit in hashed_bits], dtype=np.int8)
            return fingerprint
        
        except Exception as e:
            print(f"Error processing SMILES '{smiles}': {e}")
            return None

    def _calculate_mapc_fp(self, smiles: str) -> np.array:
        """Calculate MAP chiral fingerprints for a single SMILES string"""
        return np.array(mapc_enc(Chem.MolFromSmiles(smiles), max_radius=self.radius, n_permutations=self.fp_size, mapping=False))
        
    def _fp_method(self):
        # TODO: Add support for map4
        """
        Returns the correct function to calculate the fingerprints based on user selection
        """
        if self.fingerprint_type == 'mhfp':
            return self._calculate_mhfp_fp 
        elif self.fingerprint_type == 'mqn':
            return self._calculate_mqn_fp
        elif self.fingerprint_type == 'mapc':
            return self._calculate_mapc_fp
        elif self.fingerprint_type == 'morgan':
            return self._calculate_morgan_fp                
    
    def calculate_fingerprints(self)-> npt.NDArray:
        with Pool(processes=N_JOBS) as pool: 
            fingerprints = pool.map(self._fp_method(), self.smiles_list)

        return np.array(fingerprints)

"""
Example of usage: 
## Initiate 
fp_calculator = FingerprintCalculator(your_dataframe['smiles'], 'mqn') --> This will calculate the MQN fingerprints. For example is used when running the clustering pipeline
fp_calculator = FingerprintCalculator(smiles_list, 'mhfp', permutations=1024) --> This will calculate the MHFP with 1024 permutations (default is 512)

## Calculate fingerprints
fingerpints = fp_calculator.calculate_fingerprints() -> returns np.array witht the fp vectors of same length as smiles_list
"""
