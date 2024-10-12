import os

# =======================
# General Configuration
# =======================

# Base directory of the project
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# =======================
# File Paths
# =======================

# Input data file path
DATA_FILE_PATH = os.path.join(BASE_DIR, 'data', 'cleaned_dataset.csv')

# Output data file path
OUTPUT_FILE_PATH = os.path.join(BASE_DIR, 'data', 'processed_data.csv')

# =======================
# Loading parameters
# =======================
CHUNKSIZE = 512


# =======================
# PCA Parameters
# =======================

# Number of components for PCA
PCA_N_COMPONENTS = 3

# PCA Batch Size (if using Incremental PCA)
PCA_BATCH_SIZE = 10000

# =======================
# Logging Configuration
# =======================

# Logging levels: DEBUG, INFO, WARNING, ERROR, CRITICAL
LOGGING_LEVEL = 'INFO' #sets the threshold for the logger 

# Logging format
LOGGING_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# Log file path (optional) 
LOG_FILE_PATH = os.path.join(BASE_DIR, 'logs', 'app.log')

# =======================
# Other Configurations
# =======================


# Number of CPU cores to use for parallel processing
N_JOBS = os.cpu_count()  # Uses all available cores

# Random seed for reproducibility
RANDOM_STATE = 42
