# File to run the Clustering pipeline
import os
import glob 
import sys
import time
import logging 
import argparse
import shutil
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.decomposition import IncrementalPCA
from memory_profiler import profile
import pandas as pd

# This allows Python to recognize the src/ directory as a package and facilitates absolute imports from it.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import configurations and modules
from config import (INPUT_TMAP_PATH, OUTPUT_TMAP_PATH, CLUSTER_DATA_PATH, BASE_DIR,
                    LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH, N_JOBS, TMAP_K, TMAP_NAME)
from src.tmap_generator import TmapGenerator, ClickhouseTMAP




def format_time(seconds):
      hours, rem = divmod(seconds, 3600)
      minutes, seconds = divmod(rem, 60)
      return f"{int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process fingerprints with flexible options.")
    
    # Optional arguments to override config settings
    parser.add_argument('--l', type=str, help="'primary' to create a primary TMAP (cluster representatives) or a secondary TMAP (indicated by cluster_id = int_int_int")
    parser.add_argument('--data-file', type=str, default=INPUT_TMAP_PATH, help="Input data file path.")
    parser.add_argument('--log', type=bool, default=False, help="Saving logs to output.log file. Default False")
    parser.add_argument('--log-level', type=str, default=LOGGING_LEVEL,
                                choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Logging verbosity level.")
    parser.add_argument('--output-dir', type=str, default=OUTPUT_TMAP_PATH, help='Output directory for TMAP files (i.e. HTML and .js)')
    parser.add_argument('--tmap-name', type=str, default=TMAP_NAME, help='Name for the tmap files that will be generated.')
    parser.add_argument('--fp', type=str, default='mhfp', help='Fingerprint to use. Supported {mhfp, morgan, mapc, map4, mqn}')
    parser.add_argument('--n-jobs', type=int, default=N_JOBS, help="Number of CPU cores to use")
    return parser.parse_args()

def setup_logging(log_level, log_output):
    if log_output == True:
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, filename=LOG_FILE_PATH) # saving logs to file
    else: 
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, stream=sys.stderr)  # Direct logs to stderr (console)
    logging.info("Logging initialized")

def main() -> None:
    start = time.time()
    args = parse_arguments()

    # Set up logging based on the parsed arugments or config defaults
    setup_logging(args.log_level, args.log)

    logging.info(f"Processing file: {args.data_file}")  
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Using K={TMAP_K} (number of neighbors)")
    logging.info(f"Using {args.n_jobs} CPU cores") 

    # Check for TMAP level
    if not args.l:
        raise ValueError("No TMAP level selected. With flag --l pass 'primary' for a primary TMAP or provide a cluster_id in the format 'int_int_int' for a secondary TMAP")

    if args.l == 'primary':
        # This line of code generates a simple TMAP
        # All configuration should be done passing the arguments either with config.py file or args.parser. 
        # Generate representative cluster. 

        tmap_generator = TmapGenerator(args.data_file, fingerprint_type=args.fp,categ_cols=['cluster_id'], output_name=args.tmap_name)
        # tmap_generator.tmap_little() # -> This will generate the TMAP from the SMILES
        tmap_generator.tmap_from_vectors() # -> This method uses a KNN for generating the TMAP from the PCA coordinates

    else:
        import clickhouse_connect 

        client = clickhouse_connect.get_client(host='localhost', port=8123)
        # all_cluster_ids = client.query("SELECT DISTINCT cluster_id FROM clustered_enamine")
        # all_cluster_ids = [cluster[0] for cluster in all_cluster_ids.result_rows] # Transform into lists from list of tuples

        clickhouse_tmap = ClickhouseTMAP() 

        # max_workers = 8  # or 8, or however many CPU cores you want to use

        # # Run in parallel
        # with ProcessPoolExecutor(max_workers=max_workers) as executor:
        #     logging.info("Starting Process")
        #     # executor.map schedules each cluster_id on a separate process (up to max_workers)
        #     # This returns an iterator of whatever the method returns.
        #     results = executor.map(clickhouse_tmap.generate_tmap_from_cluster_id, all_cluster_ids, chunksize=20)
        #     print("Batch done")
        # # If your method returns something, you can iterate over it here
        # for r in results:
        #     logging.info(f"Result: {r}")

        # clickhouse_tmap.generate_tmap_from_cluster_id(args.l)
        total = 0
        # for cluster_id in all_cluster_ids:
        start = time.time()
        clickhouse_tmap.generate_tmap_from_cluster_id(args.l, client)
        end = time.time() 
        print(f"Time per TMAP: {format_time(end - start)}")
        end
if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    start_time = time.time()
    main()
    logging.info("TMAP successfully generated")

    # # Define paths
    source_pattern = os.path.join('/home/afloresep/work/chelombus/dynamic_tmap', '*.html')
    destination_dir = os.path.join('/home/afloresep/work/chelombus/maps', 'static')
    os.makedirs(destination_dir, exist_ok=True)
    logging.info(f"Source pattern for HTML files: {source_pattern}")
    logging.info(f"Destination directory: {destination_dir}")
    # Move HTML files
    for file_path in glob.glob(source_pattern):
        logging.info(f"Moving HTML file: {file_path} to {destination_dir}")
        shutil.move(file_path, destination_dir)
    # Move JavaScript files
    source_pattern_js = os.path.join('/home/afloresep/work/chelombus/dynamic_tmap', '*.js')
    for file_path in glob.glob(source_pattern_js):
        logging.info(f"Moving JS file: {file_path} to {destination_dir}")
        shutil.move(file_path, destination_dir)

    end_time = time.time()
    logging.info("All files moved successfully")
    print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")