import os
import sys
import time
import logging 
import argparse
import h5py
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.decomposition import IncrementalPCA

# Append the parent directory to sys.path to import modules from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import configurations and modules
from config import (DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, PCA_N_COMPONENTS,
                    LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH, N_JOBS, STEPS_LIST)
from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.output_generator import OutputGenerator 
from src.dimensionality_reducer import DimensionalityReducer, get_percentiles

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process fingerprints with flexible options.")
    
    # Optional arguments to override config settings
    parser.add_argument('--data-file', type=str, default=DATA_FILE_PATH, help="Input data file path.")
    parser.add_argument('--output-dir', type=str, default=OUTPUT_FILE_PATH, help="Output data directory.")
    parser.add_argument('--chunksize', type=int, default=CHUNKSIZE, help="Chunk size for loading data.")
    parser.add_argument('--pca-components', type=int, default=PCA_N_COMPONENTS, help="Number of PCA components.")
    parser.add_argument('--log-level', type=str, default=LOGGING_LEVEL,
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Logging verbosity level.")
    parser.add_argument('--n-jobs', type=int, default=N_JOBS, help="Number of CPU cores to use.")
    parser.add_argument('--resume', type=int, default=0, help="Resume from a specific chunk.")
    parser.add_argument('--load', type=int, default=1, choices=[0, 1],
                        help="Set to 1 for faster processing (more memory), or 0 for lower memory usage (slower). Default is 1.")
    parser.add_argument('--fit', type=int, default=0, choices=[0, 1], 
                        help="Set 1 to perform the fitting of coordinates within the range defined by the percentiles 0.01 and 99.99. Default is 0. Should only be used for small datasets as it increases computational time")
    parser.add_argument('--log', type=bool, default=False, help="Saving logs to output.log file. Default False")
    return parser.parse_args()

def setup_logging(log_level, log_output):
    if log_output == True:

        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, filename=LOG_FILE_PATH) # saving logs to file
    
    else: 
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, stream=sys.stderr)  # Direct logs to stderr (console)
    
    logging.info("Logging initialized")

def process_chunk(idx, chunk, data_handler, fp_calculator, output_dir):
    """
    Process a single chunk of data by calculating fingerprints and saving them.
    This function runs in parallel using ProcessPoolExecutor.
    """
    try:
        # Check if chunk already exists (HDF5 format)
        fp_chunk_path = os.path.join(output_dir, f'fp_chunks/fingerprints_chunk_{idx}.h5')
        if os.path.exists(fp_chunk_path):
            # logging.info(f'Chunk {idx} already processed, skipping.')
            return

        # Extract smiles and features from chunk
        smiles_list, features = data_handler.extract_smiles_and_features(chunk)

        # Calculate fingerprints
        fingerprints = fp_calculator.calculate_fingerprints(smiles_list)

        # Ensure output directories exist
        os.makedirs(os.path.join(output_dir, 'fp_chunks'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'features_chunks'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'output'), exist_ok=True)

        # Save fingerprints in HDF5 format
        with h5py.File(fp_chunk_path, 'w') as h5f:
            h5f.create_dataset('fingerprints', data=fingerprints)

         # Save smiles and features in HDF5 format
        with h5py.File(os.path.join(output_dir, f'features_chunks/smiles_features_chunk_{idx}.h5'), 'w') as h5f:
            h5f.create_dataset('smiles_list', data=np.array(smiles_list, dtype='S'))  # Store strings as bytes
            # h5f.create_dataset('features', data= np.array(features, dtype='S')) 

        del fingerprints, smiles_list, features 
        
    except Exception as e:
        logging.error(f"Error processing chunk {idx}: {e}", exc_info=True)

def chunk_generator(data_chunks):
    """
    Generator to yield chunks one by one to minimize memory usage.
    """
    for chunk in data_chunks:
        yield chunk

def main():
    args = parse_arguments()

    # Set up logging based on the parsed arguments or config defaults
    setup_logging(args.log_level, args.log)
    
    logging.info(f"Processing file: {args.data_file}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Chunk size: {args.chunksize}")
    logging.info(f"Number of PCA components: {args.pca_components}")
    logging.info(f"Using {args.n_jobs} CPU cores")
    logging.info(f"Load mode: {'High memory usage (fast)' if args.load == 1 else 'Low memory usage (slow)'}")
    
    assert args.pca_components == len(STEPS_LIST), "STEPS_LIST should be same lenght as number of PCA_COMPONENTS" 

    # Initialize classes
    data_handler = DataHandler(args.data_file, args.chunksize)
    output_gen = OutputGenerator()
    fp_calculator = FingerprintCalculator()
    dim_reducer = DimensionalityReducer()
       
    # Load data in chunks
    data_chunks, total_chunks = data_handler.load_data()

    # Process chunks with tqdm progress bar
    start = time.time()

    # Process chunks in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.n_jobs) as executor:
        futures = []
        for idx, chunk in enumerate(tqdm(data_chunks, total=total_chunks, desc="Loading chunks and calculating fingerprints")):
            # Resume from a specific chunk if needed
            if idx < args.resume:
                continue
            futures.append(executor.submit(process_chunk, idx, chunk, data_handler, fp_calculator, args.output_dir))

        # Wait for all futures to complete with error handling
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing chunks"):
            try:
                future.result()  # Raises any exceptions from the workers
            except Exception as e:
                logging.error(f"Error in parallel processing: {e}", exc_info=True)

    end = time.time()
    logging.info(f"Preprocessing of data took: {(end - start)/60:.2f} minutes")

    # Partial fit using Incremental PCA
    ipca = IncrementalPCA(n_components=args.pca_components)  # Dimensions to reduce to

    for idx in tqdm(range(args.resume, total_chunks), desc="Loading Fingerprints and fitting"):
        try:
            fp_chunk_path = os.path.join(args.output_dir, f'fp_chunks/fingerprints_chunk_{idx}.h5')
            with h5py.File(fp_chunk_path, 'r') as h5f:
                fingerprints = h5f['fingerprints'][:]

            ipca.partial_fit(fingerprints)

            del fingerprints

        except Exception as e:
            logging.error(f"Error during PCA fitting for chunk {idx}: {e}", exc_info=True)

    # Transform Data and Save results
    logging.info("Performing Dimensionality Reduction...")

    for idx in tqdm(range(args.resume, total_chunks), desc='Transforming Data'):
        try:
            # Load fingerprint
            fp_chunk_path = os.path.join(args.output_dir, f'fp_chunks/fingerprints_chunk_{idx}.h5')
            with h5py.File(fp_chunk_path, 'r') as h5f:
                fingerprints = h5f['fingerprints'][:]

            # Load smiles and features
            feat_chunk_path = os.path.join(args.output_dir, f'features_chunks/smiles_features_chunk_{idx}.h5')
            with h5py.File(feat_chunk_path, 'r') as h5f:
                smiles_list = h5f['smiles_list'][:].astype(str)  # Convert from bytes to strings
                # features = h5f['features'][:].astype(str)
                features = []

            coordinates = ipca.transform(fingerprints)
              
            # Output coordinates before clipping with Percentiles.
            output_gen.save_batch(idx, coordinates, smiles_list, features, args.output_dir)

            if args.fit == 1:
               digest_methods = dim_reducer.digest_generator(args.pca_components)

               for i in range(len(digest_methods)):
                    digest_methods[i].batch_update(coordinates[:, i])

            # Free memory space
            del fingerprints, coordinates, smiles_list, features 
            
            # Delete intermediate files after loading them 
            os.remove(fp_chunk_path)
            os.remove(feat_chunk_path)

        except Exception as e:
            logging.error(f"Error during data transformation for chunk {idx}: {e}", exc_info=True)


    # Get Percentiles 
    if args.fit == 1:
        try:
            percentiles = get_percentiles(digest_methods, STEPS_LIST)
            print(percentiles)

        except Exception as e:
            logging.error(f"Error calculating percentiles: {e}", exc_info=True)
            # Use default percentiles if calculation fails?

         # Map PCA coordinates 
        output_gen.steps = [32, 32, 32]  # Number of steps to be taken on each dimension
        start = time.time()

        logging.info('Fitting coordinates to cube')
        try:
             for output_file in os.listdir(os.path.join(args.output_dir, 'output')):
                 output_gen.fit_coord_multidimensional(output_file, percentiles)

        except Exception as e:
             logging.error(f"Error fitting coordinates to cube: {e}", exc_info=True)

        end = time.time()
        logging.info(f'Total time fitting coordinates: {(end - start)/60:.2f} minutes')

        
if __name__ == '__main__':
    start_time = time.time() 
    main() 
    end_time = time.time()
    logging.info(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")
