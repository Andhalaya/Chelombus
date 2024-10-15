import os
import time
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
from config import DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, N_JOBS
from tqdm import tqdm  # Import tqdm
from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.output_generator import OutputGenerator
import time

# Define the chunk sizes you want to test
data_handler = DataHandler()
fp_calculator = FingerprintCalculator()
chunk_sizes = range(5000, 6650,150)  # Example sizes
timing_results = {}

for chunk_size in chunk_sizes:
    start_time = time.time()

    data_chunks = data_handler.load_data('/home/afloresep/work/chelombus/data/10M_ZINC_id_Sim_mqn.csv', chunk_size)

    for idx, chunk in enumerate(data_chunks):
        smiles_list, features = data_handler.extract_smiles_and_features(chunk)
        fingerprints = fp_calculator.calculate_fingerprints(smiles_list)

    end_time = time.time()
    elapsed_time = end_time - start_time
    timing_results[chunk_size] = elapsed_time/60
    print(f"Processing with chunk size {chunk_size} took {elapsed_time/60:.2f} minutes.")

print("Timing results:", timing_results)
