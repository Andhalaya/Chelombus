import os
import time
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
from config import DATA_FILE_PATH, CHUNKSIZE, PCA_N_COMPONENTS
from tqdm import tqdm 
from sklearn.decomposition import IncrementalPCA

from src.data_handler import DataHandler, get_total_chunks
from src.fingerprint_calculator import FingerprintCalculator
from src.output_generator import OutputGenerator
import numpy as np 
from tdigest import TDigest

def main():
    # Initialize classes
    data_handler = DataHandler(DATA_FILE_PATH, 36050)
    output_gen = OutputGenerator()
    fp_calculator = FingerprintCalculator()

    # # TODO: Dynamic method for calculating the percentiles for every PCA dimension 
    # digests = {}
    # for i in range(5):
    #     digests[f'{i}_digest'] = TDigest()

    x_digest = TDigest()
    y_digest = TDigest()
    z_digest = TDigest()
    
    # Load data in chunks
    data_chunks, total_chunks = data_handler.load_data()

    # Process chunks with tqdm progress bar

    start = time.time()

    for idx, chunk in enumerate(tqdm(data_chunks, total=total_chunks, desc="Loading chunks and calculating fingerprints")):

        # Check if chunk already exists
        if os.path.exists(f'data/10M/fp_chunks/fingerprints_chunk_{idx}.pkl'):
            continue
        
        # Extract smiles and features from chunk
        smiles_list, features = data_handler.extract_smiles_and_features(chunk)

        # Calculate fingerprints with progress bar
        fingerprints = fp_calculator.calculate_fingerprints(smiles_list)

        # Save  fingerprints
        with open(f'data/1M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((fingerprints), f)

        # Save rest of data
        with open(f'data/1M/features_chunks/smiles_features_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((smiles_list, features), f)

        del smiles_list, features, fingerprints # Free space

    end = time.time()

    print(f"Preprocessing of data took: {(end-start)/60} mninutes")


    # Partial fit using iPCA
    ipca = IncrementalPCA(n_components = PCA_N_COMPONENTS) # Dimensions to reduce to
    
    for idx in tqdm(range(total_chunks), desc="Loading Fingerprints and fitting "):
        with open(f'data/1M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        ipca.partial_fit(fingerprints)
        del fingerprints

    # Transform Data and Save results
    print("Performing Dimensionality Reduction...")
    for idx in tqdm(range(total_chunks), desc='Transforming Data'):
        # Load fingerprint
        with open(f'data/1M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        coordinates = ipca.transform(fingerprints)

        x_digest.batch_update(coordinates[:,0])
        y_digest.batch_update(coordinates[:,1])
        z_digest.batch_update(coordinates[:,2])

        with open(f'data/1M/features_chunks/smiles_features_chunk_{idx}.pkl', 'rb') as f:
            smiles_list, features = pickle.load(f)

        output_gen.save_batch(idx, coordinates, smiles_list, features)

        # Clean
        del fingerprints, coordinates, smiles_list, features
    
    print(f'X: {x_digest.percentile(0.01), x_digest.percentile(99.99)}')
    print(f'Y: {y_digest.percentile(0.01), x_digest.percentile(99.99)}')
    print(f'Z: {z_digest.percentile(0.01), x_digest.percentile(99.99)}')

if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print(f"Total time: {(end - start)/60} minutes")
