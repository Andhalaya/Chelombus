# File to run the Clustering pipeline

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

# This allows Python to recognize the src/ directory as a package and facilitates absolute imports from it.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import configurations and modules
from config import (DATA_FILE_PATH, INPUT_TMAP_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, OUTPUT_TMAP_PATH, PCA_N_COMPONENTS,
                    LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH, N_JOBS, STEPS_LIST)
from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.tmap_generator import TmapConstructor, TmapGenerator

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process fingerprints with flexible options.")
    
    # Optional arguments to override config settings
    parser.add_argument('--data-file', type=str, default=INPUT_TMAP_PATH, help="Input data file path.")
    parser.add_argument('--log', type=bool, default=False, help="Saving logs to output.log file. Default False")
    parser.add_argument('--log-level', type=str, default=LOGGING_LEVEL,
                                choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Logging verbosity level.")
    parser.add_argument('--output-dir', type=str, default=OUTPUT_TMAP_PATH, help='Output directory for TMAP files (i.e. HTML and .js)')
    parser.add_argument('--n-jobs', type=int, default=N_JOBS, help="Number of CPU cores to use")
    return parser.parse_args()

def setup_logging(log_level, log_output):
    if log_output == True:
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, filename=LOG_FILE_PATH) # saving logs to file
    else: 
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, stream=sys.stderr)  # Direct logs to stderr (console)
    logging.info("Logging initialized")

def main() -> None:
    args = parse_arguments()

    # Set up logging based on the parsed arugments or config defaults
    setup_logging(args.log_level, args.log)
    main_dataframe = pd.read_csv(INPUT_TMAP_PATH) 
    fp_calculator = FingerprintCalculator(main_dataframe['smiles'], fingerprint_type='morgan',fp_size=1024)
    tmap_constructor = TmapConstructor(main_dataframe) 
    tmap_generator = TmapGenerator(main_dataframe, 'morgan', categ_cols=['cluster_id'])
    logging.info(f"Processing file: {args.data_file}")  
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Using {args.n_jobs} CPU cores") 
    
    logging.info(f"Calculating fingerprints...") 
    fingerprints = fp_calculator.calculate_fingerprints()

    logging.info("Constructing LSH Forest")
    tmap_generator.construct_lsh_forest()

    logging.info("Plotting...")
    properties = tmap_constructor.mol_properties_from_df()
    print(len(properties))
    # tmap_generator.plot_faerun()


if __name__ == '__main__':
    start_time = time.time() 
    main() 
    end_time = time.time()
    logging.info(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")
