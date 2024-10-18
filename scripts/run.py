import os
import time
import sys
import pickle
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from config import DATA_FILE_PATH, CHUNKSIZE, PCA_N_COMPONENTS
from tqdm import tqdm 
from sklearn.decomposition import IncrementalPCA

from src.data_handler import DataHandler, get_total_chunks
from src.fingerprint_calculator import FingerprintCalculator
from src.output_generator import OutputGenerator, get_percentiles
import numpy as np 
import inspect
from tdigest import TDigest 

def main():
    # Initialize classes
    data_handler = DataHandler(DATA_FILE_PATH, 36050)
    output_gen = OutputGenerator()
    fp_calculator = FingerprintCalculator()
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
        with open(f'data/10M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((fingerprints), f)

        # Save rest of data
        with open(f'data/10M/features_chunks/smiles_features_chunk_{idx}.pkl', 'wb') as f:
            pickle.dump((smiles_list, features), f)

        del smiles_list, features, fingerprints # Free space

    end = time.time()

    print(f"Preprocessing of data took: {(end-start)/60} minutes")


    # Partial fit using iPCA
    ipca = IncrementalPCA(n_components = PCA_N_COMPONENTS) # Dimensions to reduce to
    
    for idx in tqdm(range(total_chunks), desc="Loading Fingerprints and fitting "):
        with open(f'data/10M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        ipca.partial_fit(fingerprints)
        del fingerprints

    # Transform Data and Save results
    print("Performing Dimensionality Reduction...")

    for idx in tqdm(range(total_chunks), desc='Transforming Data'):
        # Load fingerprint
        with open(f'data/10M/fp_chunks/fingerprints_chunk_{idx}.pkl', 'rb') as f:
            fingerprints = pickle.load(f)

        # Load smiles and features
        with open(f'data/10M/features_chunks/smiles_features_chunk_{idx}.pkl', 'rb') as f:
            smiles_list, features = pickle.load(f)

        # Get coordinates in np.array(n_smiles, n_pca_dim)
        coordinates = ipca.transform(fingerprints)
        
        # Update digest for every batch
        x_digest.batch_update(coordinates[:,0])
        y_digest.batch_update(coordinates[:,1])
        z_digest.batch_update(coordinates[:,2])

        # Output coordiantes ~before clip with Percentiles. 
        output_gen.save_batch(idx, coordinates, smiles_list, features)
    
        del fingerprints, coordinates, smiles_list, features
    
    percentiles = get_percentiles(x_digest,  y_digest, z_digest)

    # percentiles_10M = [(-23.69173900049878, 29.46851433296263), (-16.738348923448807, 23.54345794491185), (-10.819449416425204, 12.994996725530592)] 
    
    # Mapp PCA coordinates to the 100x100x100 dimensional cube
    # mapped_coordinates = output_gen.map_to_grid(coordinates, percentiles)

    # Output 
    

    # Clean


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print(f"Total time: {(end - start)/60} minutes")