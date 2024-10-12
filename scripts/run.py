import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
from config import DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, N_JOBS
from tqdm import tqdm  # Import tqdm

from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.dimensionality_reducer import DimensionalityReducer
from src.output_generator import OutputGenerator

def main():
    # Initialize classes
    data_handler = DataHandler()
    fp_calculator = FingerprintCalculator()
    reducer = DimensionalityReducer()
    output_gen = OutputGenerator()

    # Calculate total number of chunks (optional, for tqdm)
    total_chunks = data_handler.get_total_chunks(DATA_FILE_PATH, CHUNKSIZE)

    # Load data in chunks
    data_chunks = data_handler.load_data(DATA_FILE_PATH, CHUNKSIZE)

    # Process chunks with tqdm progress bar
    # Partial fit using iPCA 
    num_chunks = 0
    for idx, chunk in enumerate(tqdm(data_chunks, total=total_chunks, desc="Processing Chunks")):
        num_chunks += 1

        # Check if chunk already exists
        if os.path.exists(f'data/fingerprints_chunk_{idx}.pkl'):
            pass
        else:
            smiles_list, features = data_handler.extract_smiles_and_features(chunk)

            # Calculate fingerprints with progress bar
            fingerprints = fp_calculator.calculate_fingerprints(smiles_list)
            
            # Save fingerprints
            with open(f'data/fingerprints_chunk_{idx}.pkl', 'wb') as f:
                pickle.dump((fingerprints, smiles_list, features), f)

    # Load all fingerprints with tqdm progress bar
    fingerprints, smiles_list, oh_features = [], [], []
    for idx in tqdm(range(num_chunks), desc="Loading Fingerprints"):
        with open(f'data/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fps_chunk, smiles_chunk, features_chunk = pickle.load(f)
            fingerprints.extend(fps_chunk)
            smiles_list.extend(smiles_chunk)
            oh_features.extend(features_chunk)

    # data_dict = { 'smiles': smiles_list, 'fp': fingerprints, 'oh_features':oh_features} 
    # # Perform Dimensionality Reduction
    # print("Performing Dimensionality Reduction...")
    # coordinates = reducer.fit_transform(fingerprints)

    # # Compile and export data
    # output_gen.compile_output(coordinates, smiles_list, features)
    # output_gen.export_data(OUTPUT_FILE_PATH)

if __name__ == '__main__':
    main()
