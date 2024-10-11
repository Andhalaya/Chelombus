from multiprocessing import Pool
from config import N_JOBS

class FingerprintCalculator:
    def calculate_fingerprints(self, smiles_list):
        with Pool(processes=N_JOBS) as pool:
            fingerprints = pool.map(self._compute_fingerprint, smiles_list)
        return fingerprints
    
    # handle exception for invalid SMILES string
    # Process data in chunks or using generators? 