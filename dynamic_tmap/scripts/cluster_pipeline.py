# File to run the Clustering pipeline

import os
import sys
import time
import logging 
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.decomposition import IncrementalPCA
from memory_profiler import profile
import gc
import pandas as pd

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
    parser.add_argument('--log', type=bool, default=False, help="Saving logs to output.log file. Default False")
    return parser.parse_args()

def setup_logging(log_level, log_output):
    if log_output == True:

        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, filename=LOG_FILE_PATH) # saving logs to file
    
    else: 
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, stream=sys.stderr)  # Direct logs to stderr (console)
    
    logging.info("Logging initialized")


def main() -> None:
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
    dim_reducer = DimensionalityReducer()
       
    # Get correct generator function for loading chunks and number of chunks for tqmd 
    data_chunk_loader, total_chunks = data_handler.load_data()

    #TODO: Add GPU fingerprint calculation support
    # Process each chunk: Data chunk loader will yield a chunk and process_chunk will calculate its fingerprints & save them to parquet files
    start = time.time()
    for idx, chunk in enumerate(tqdm(data_chunk_loader, total= total_chunks, desc=f"Loading chunk and calculating its fingerprints")):
        data_handler.process_chunk(idx, chunk, args.output_dir)
        del idx, chunk
    end = time.time()    
    logging.info(f"Preprocessing of data took: {(end - start)/60:.2f} minutes")
 
    # Incremental PCA fitting
    ipca = IncrementalPCA(n_components=args.pca_components)  # Dimensions to reduce to
    
    for idx in tqdm(range(args.resume, total_chunks), desc="Loading Fingerprints and iPCA partial fitting"):
        try:
            fp_chunk_path = os.path.join(args.output_dir, f'batch_parquet/fingerprints_chunk_{idx}.parquet')
            df_fingerprints = pd.read_parquet(fp_chunk_path, engine="pyarrow")
            
            data = df_fingerprints.drop(columns=['smiles']).values
            ipca.partial_fit(data)

            del df_fingerprints, data
            gc.collect()

        except Exception as e:
            logging.error(f"Error during PCA fitting for chunk {idx}: {e}", exc_info=True)

    # Transform Data and Save results
    logging.info("Performing Dimensionality Reduction...")

    for idx in tqdm(range(args.resume, total_chunks), desc='iPCA transform and saving results'):
        try:
            # Load fingerprint
            fp_chunk_path = os.path.join(args.output_dir, f'batch_parquet/fingerprints_chunk_{idx}.parquet')
            df_fingerprints  =  pd.read_parquet(fp_chunk_path)
            features = [] # Is this necessary?  
            coordinates = ipca.transform(df_fingerprints.drop(columns=['smiles']).values)  # -> np.array shape (chunk_size, n_pca_comp)

            # Output coordinates into a parquet file.
            output_gen.batch_to_multiple_parquet(idx, coordinates,df_fingerprints['smiles'].to_list(), features, args.output_dir)

            # Free memory space
            del df_fingerprints, coordinates, features 
            # os.remove(fp_chunk_path) 

        except Exception as e:
            logging.error(f"Error during data transformation for chunk {idx}: {e}", exc_info=True)
        
if __name__ == '__main__':
    start_time = time.time()
    # Run the main function
    main()
    # Manually invoke garbage collection
    gc.collect()
    # Clear all modules (optional but useful for memory cleanup)
    sys.modules.clear()
    # Clear any cached variables or objects
    del main
    end_time = time.time()
    logging.info(f"Total execution time: {int((end_time - start_time) // 3600)} hours, {int(((end_time - start_time) % 3600) // 60)} minutes, and {int((end_time - start_time) % 60)} seconds")