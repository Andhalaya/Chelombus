import os
import time
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
from config import DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, N_JOBS
from tqdm import tqdm  # Import tqdm
from sklearn.decomposition import IncrementalPCA

from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.dimensionality_reducer import DimensionalityReducer
from src.output_generator import OutputGenerator

def main():
    # Initialize classes
    data_handler = DataHandler()
    fp_calculator = FingerprintCalculator()
    # reducer = DimensionalityReducer()
    output_gen = OutputGenerator()

    # Calculate total number of chunks (optional, for tqdm)
    start = time.time() 
    total_chunks = data_handler.get_total_chunks(DATA_FILE_PATH, CHUNKSIZE)
    end = time.time()

    print(f"total chunks calculated in: {end-start} seconds")

    # Load data in chunks
    start = time.time()
    data_chunks = data_handler.load_data(DATA_FILE_PATH, CHUNKSIZE)

    # Process chunks with tqdm progress bar
    num_chunks = 0
    for idx, chunk in enumerate(tqdm(data_chunks, total=total_chunks, desc="Processing Chunks")):
        num_chunks += 1

        # Check if chunk already exists
        if os.path.exists(f'data/fingerprints_chunk_{idx}.pkl'):
            continue

        # Load smiles and features
        smiles_list, features = data_handler.extract_smiles_and_features(chunk)

        # Calculate fingerprints with progress bar
        fingerprints = fp_calculator.calculate_fingerprints(smiles_list)

        # Save fingerprints
        with open(f'data/fp_chunks/fingerprints_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((fingerprints), f)

        # Save rest of data
        with open(f'data/features_chunks/smiles_features_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((smiles_list, features), f)

        del smiles_list, features, fingerprints # Free space

    end = time.time()
    print(f"All fingerprints calculated in: {end-start} seconds")


    # Load all fingerprints with tqdm progress bar
    # Partial fit using iPCA
    ipca = IncrementalPCA(n_components = 3) # Dimensions to reduce to
    for idx in tqdm(range(num_chunks), desc="Loading Fingerprints and fitting "):
        with open(f'data/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        ipca.partial_fit(fingerprints)
        del fingerprints

    # Transform Data and Save results
    print("Performing Dimensionality Reduction...")

    for idx in tqdm(range(total_chunks), desc='Transforming Data'):
        # Load fingerprint
        with open(f'data/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        coordinates = ipca.fit_transform(fingerprints)

        with open(f'data/smiles_features_chunk_{idx}.pkl', 'rb') as f:
            smiles_list, features = pickle.load(f)
        
        output_gen.save_batch(idx, coordinates, smiles_list, features)

        # Clean
        del fingerprints, coordinates, smiles_list, features
    
if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print(f"Total time: {(end - start)/60} minutes")
