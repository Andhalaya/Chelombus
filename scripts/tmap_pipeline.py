# File to run the TMAP pipeline

# File to run the Clustering pipeline

import os
import sys
import time
import logging 
import argparse
import numpy as np
from tqdm import tqdm
import pandas as pd

# Append the parent directory to sys.path to import modules from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.fingerprint_calculator import FingerprintCalculator
from src.layout_computer import LayoutComputer
from src.tmap_generator import TmapConstructor

# Import configurations and modules
from config import (DATA_FILE_PATH, OUTPUT_FILE_PATH, CHUNKSIZE, PCA_N_COMPONENTS,
                    LOGGING_LEVEL, LOGGING_FORMAT, LOG_FILE_PATH, N_JOBS, STEPS_LIST)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process fingerprints with flexible options.")
    
    # Optional arguments to override config settings
    parser.add_argument('--data-file', type=str, default=DATA_FILE_PATH, help="Input data file path.")
    parser.add_argument('--log', type=bool, default=False, help="Saving logs to output.log file. Default False")
    return parser.parse_args()

def setup_logging(log_level, log_output):
    if log_output == True:
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, filename=LOG_FILE_PATH) # saving logs to file
    else: 
        logging.basicConfig(level=getattr(logging, log_level.upper()), format=LOGGING_FORMAT, stream=sys.stderr)  # Direct logs to stderr (console)
    logging.info("Logging initialized")

def main() -> None:
    pass