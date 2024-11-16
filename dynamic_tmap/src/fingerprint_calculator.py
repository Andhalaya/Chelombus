from multiprocessing import Pool
from functools import partial
from typing import Optional
from config import N_JOBS
from rdkit import Chem
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import numpy as np
import tmap as tm
import hashlib
import cudf
import cupy as cp
from concurrent.futures import ThreadPoolExecutor
from mapchiral.mapchiral import encode as mapc_enc
import numpy.typing as npt

# Define the fingerprint functions at the module level

def calculate_mhfp_fp_gpu(smiles: str, permutations: int) -> cp.ndarray:
    """Calculate MHFP fingerprint for a single SMILES string using GPU."""
    try:
        encoder = MHFPEncoder(permutations)
        return cp.array(encoder.encode(smiles))
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return cp.zeros(permutations)

def calculate_mqn_fp_gpu(smiles: str) -> cp.ndarray:
    """Calculate MQN fingerprint for a single SMILES string using GPU."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return cp.zeros(42)  # MQN fingerprint has 42 features
        fingerprint = rdMolDescriptors.MQNs_(mol)
        return cp.array(fingerprint)
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return cp.zeros(42)

def calculate_morgan_fp_gpu(smiles: str, radius: int, fp_size: int) -> cp.ndarray:
    """Calculate Morgan fingerprint for a single SMILES string using GPU."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string.")
        morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=fp_size)
        return cp.array(morgan_fp)
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return cp.zeros(fp_size)

def calculate_mapc_fp_gpu(smiles: str, radius: int, fp_size: int) -> cp.ndarray:
    """Calculate MAP Chiral fingerprint for a single SMILES string using GPU."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string.")
        fingerprint = mapc_enc(mol, max_radius=radius, n_permutations=fp_size, mapping=False)
        return cp.array(fingerprint)
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return cp.zeros(fp_size)

class FingerprintCalculatorGPU:
    def __init__(
        self,
        smiles_list: list,
        fingerprint_type: str,
        permutations: Optional[int] = 512,
        fp_size: Optional[int] = 1024,
        radius: Optional[int] = 2,
        chunk_size: Optional[int] = 100_000
    ):
        self.smiles_list = smiles_list
        self.fingerprint_type = fingerprint_type.lower()
        self.permutations = permutations
        self.fp_size = fp_size
        self.radius = radius
        self.chunk_size = chunk_size

        # Map fingerprint types to functions
        self.fingerprint_function_map = {
            'mhfp': calculate_mhfp_fp_gpu,
            'mqn': calculate_mqn_fp_gpu,
            'mapc': calculate_mapc_fp_gpu,
            'morgan': calculate_morgan_fp_gpu,
        }

    def parse_smiles_in_parallel(self, max_workers=16):
        """Parse SMILES strings using ThreadPoolExecutor for faster CPU parsing."""
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            mols = list(executor.map(Chem.MolFromSmiles, self.smiles_list))
        return mols

    def calculate_fingerprints_gpu(self) -> cp.ndarray:
        """Calculate fingerprints using GPU."""
        func_to_apply = self.fingerprint_function_map.get(self.fingerprint_type)

        if func_to_apply is None:
            raise ValueError(f"Unsupported fingerprint type: {self.fingerprint_type}")

        # Prepare the function with necessary parameters
        if self.fingerprint_type == 'mhfp':
            func_to_apply = partial(func_to_apply, permutations=self.permutations)
        elif self.fingerprint_type == 'morgan':
            func_to_apply = partial(func_to_apply, radius=self.radius, fp_size=self.fp_size)
        elif self.fingerprint_type == 'mapc':
            func_to_apply = partial(func_to_apply, radius=self.radius, fp_size=self.fp_size)
        # MQN doesn't need this
        
        fingerprints = []

        # Process in chunks to fit into GPU memory
        for i in range(0, len(self.smiles_list), self.chunk_size):
            chunk = self.smiles_list[i:i + self.chunk_size]
            mols = self.parse_smiles_in_parallel()
            chunk_fps = [func_to_apply(Chem.MolToSmiles(mol)) for mol in mols if mol is not None]
            chunk_fps_gpu = cp.asarray(chunk_fps)
            fingerprints.append(chunk_fps_gpu)

        # Concatenate all fingerprints into a single CuPy array
        all_fingerprints_gpu = cp.concatenate(fingerprints)
        return all_fingerprints_gpu

# Example usage
if __name__ == "__main__":
    smiles_list = ["CCO", "CCN", "CCC", "CCCl", "CCBr"]
    calculator = FingerprintCalculatorGPU(smiles_list, fingerprint_type='mqn')
    fingerprints_gpu = calculator.calculate_fingerprints_gpu()
    print(f"Fingerprints shape: {fingerprints_gpu.shape}")
