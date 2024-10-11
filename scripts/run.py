import os
import pickle
from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.dimensionality_reducer import DimensionalityReducer
from src.output_generator import OutputGenerator
from config import DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE


def main():
    # Initialize classes
    data_handler = DataHandler()
    fp_calculator = FingerprintCalculator()
    reducer = DimensionalityReducer()
    output_gen = OutputGenerator()

    # Load in chunks
    data_chunks = data_handler.load_data_in_chunks(DATA_FILE_PATH, CHUNKSIZE)
    num_chunks = 0

    # Process chunks
    for idx, chunk in enumerate(data_chunks): 
        num_chunks += 1
        # Check if chunk already exists
        if os.path.exists(f'data/fingerprints_chunk_{idx}.pkl'):
            print(f'Chunk {idx} already processed. Skipping.')

        else: 
            print(f'Processing chunk {idx}')
            smiles_list = data_handler.extract_smiles(chunk)
            features = data_handler.extract_features(chunk)

            fingerprints = fp_calculator.calculate_fp(smiles_list)
            
            # Save fingerprints
            with open(f'data/fingerprints_chunk_{idx}.pkl', 'wb') as f:
                pickle.dump((fingerprints, smiles_list, features), f)


    # Load all fingerprints
    fingerprints, smiles_list, features = [], [], []
    for idx in range(num_chunks):
        with open(f'data/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fps_chunk, smiles_chunk, features_chunk = pickle.load(f)
            fingerprints.extend(fps_chunk)
            smiles_list.extend(smiles_chunk)
            features.extend(features_chunk)


    # Dimensionality Reduction
    coordinates = reducer.fit_transform(fingerprints)

    # Compile and export data
    output_gen.compile_output(coordinates, smiles_list, features)
    output_gen.export_data('data/processed_data.csv')

if __name__ == '__main__':
    main()
