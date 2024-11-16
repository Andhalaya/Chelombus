import time
import random
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import cupy as cp
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from functools import partial

# Function for generating random SMILES strings (simplified version)
def generate_random_smiles(num_smiles):
    atoms = ["C", "N", "O", "F", "Cl", "Br"]
    return ["".join(random.choices(atoms, k=5)) for _ in range(num_smiles)]

# CPU-based MQN fingerprint calculation
def calculate_mqn_fp(smiles: str) -> np.ndarray:
    try:
        fingerprint = rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smiles))
        return np.array(fingerprint)
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return np.zeros(42)

# GPU-based MQN fingerprint calculation
def calculate_mqn_fp_gpu(smiles: str) -> cp.ndarray:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return cp.zeros(42)
        fingerprint = rdMolDescriptors.MQNs_(mol)
        return cp.array(fingerprint)
    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        return cp.zeros(42)

# CPU version using multiprocessing
def calculate_fingerprints_cpu(smiles_list, n_jobs=16):
    with Pool(processes=n_jobs) as pool:
        fingerprints = pool.map(calculate_mqn_fp, smiles_list)
    return np.array(fingerprints)

# GPU version using CuPy
def calculate_fingerprints_gpu(smiles_list):
    with ThreadPoolExecutor(max_workers=16) as executor:
        mols = list(executor.map(Chem.MolFromSmiles, smiles_list))
    fingerprints = [calculate_mqn_fp_gpu(Chem.MolToSmiles(mol)) for mol in mols if mol is not None]
    return cp.array(fingerprints)

# Benchmarking function
def run_benchmark(num_smiles):
    # Generate random SMILES strings
    smiles_list = generate_random_smiles(num_smiles)

    # CPU Benchmark
    print("Running CPU benchmark...")
    start_time = time.time()
    cpu_fingerprints = calculate_fingerprints_cpu(smiles_list)
    cpu_duration = time.time() - start_time
    print(f"CPU time: {cpu_duration:.2f} seconds")

    # GPU Benchmark
    print("Running GPU benchmark...")
    start_time = time.time()
    gpu_fingerprints = calculate_fingerprints_gpu(smiles_list)
    gpu_duration = time.time() - start_time
    print(f"GPU time: {gpu_duration:.2f} seconds")

    # Check for correctness
    cpu_result = cpu_fingerprints[0]
    gpu_result = cp.asnumpy(gpu_fingerprints[0])
    if np.allclose(cpu_result, gpu_result):
        print("Results match!")
    else:
        print("Results do not match!")

    # Print speedup
    speedup = cpu_duration / gpu_duration
    print(f"Speedup: {speedup:.2f}x")

# Run the benchmark with 1 million SMILES strings
if __name__ == "__main__":
    run_benchmark(1_000_000)
