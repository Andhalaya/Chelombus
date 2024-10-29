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
    parser.add_argument('--load', type=int, default=0, choices=[0, 1],
                        help="Set to 1 for faster processing (fitting entire data into memory), or 0 for lower memory usage (one chunk at a timer). Default is 1.")
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
    fp_calculator = FingerprintCalculator()
    dim_reducer = DimensionalityReducer()
       
    # Load data in chunks
    data_chunks, total_chunks = data_handler.load_data()

    # Process chunks with tqdm progress bar
    start = time.time()

    # Process chunks in parallel using ProcessPoolExecutor
    if args.load == 1:
        """
        This will process all chunks in parallel. It ensures fast fingerprint calculation at the cost of loading all data into memory
        Use only if the whole dataset fits into memory. Otherwise, use sequential chunk calculation
        """
        with ProcessPoolExecutor(max_workers=args.n_jobs) as executor:
            futures = []
            for idx, chunk in enumerate(tqdm(data_chunks, total=total_chunks, desc="Loading chunks")):
                # Resume from a specific chunk if needed
                if idx < args.resume:
                    continue
                futures.append(executor.submit(data_handler.process_chunk, idx, chunk, data_handler, fp_calculator, args.output_dir))
                del idx, chunk
            # Wait for all futures to complete with error handling

            for future in tqdm(as_completed(futures), total=len(futures), desc="Calculating fingerprints for all chunks"):
                try:
                    future.result()  # Raises any exceptions from the workers
                except Exception as e:
                    logging.error(f"Error in parallel processing: {e}", exc_info=True)

        # Ensure all process have completed before moving on with Loading Fingerprints and iPCA partial fitting
        executor.shutdown(wait=True)

    if args.load == 0:
       """
       This will instead process each chunk at a time to ensure that the entire dataset is not loaded into memory
       You can use higher chunksize in this method
       """
       for idx, chunk in enumerate(tqdm(data_chunks, total= total_chunks, desc=f"Loading chunk and calculating its fingerprints")):
           data_handler.process_chunk(idx, chunk, fp_calculator, args.output_dir)
           del idx, chunk

    end = time.time()    
    logging.info(f"Preprocessing of data took: {(end - start)/60:.2f} minutes")
 
    ipca = IncrementalPCA(n_components=args.pca_components)  # Dimensions to reduce to

    # Incremental PCA fitting
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
            features = []
            coordinates = ipca.transform(df_fingerprints.drop(columns=['smiles']).values)  # -> np.array shape (chunk_size, n_pca_comp)

            # Output coordinates into a parquet file.
            output_gen.batch_to_one_parquet(coordinates,df_fingerprints['smiles'].to_list(), features, args.output_dir)

            if args.fit == 1:
               digest_methods = dim_reducer.digest_generator(args.pca_components)

               for i in range(len(digest_methods)):
                    digest_methods[i].batch_update(coordinates[:, i])

            # Free memory space
            del df_fingerprints, coordinates, features 
            os.remove(fp_chunk_path) 

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
