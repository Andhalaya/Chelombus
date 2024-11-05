import os

# =======================
# CLUSTER Configuration
# =======================

# Base directory of the project
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# =======================
# File Paths
# =======================

# Input data file path
DATA_FILE_PATH = os.path.join(BASE_DIR, 'data', '230M.cxsmiles')

# Output data file path
OUTPUT_FILE_PATH = os.path.join(BASE_DIR, 'data/230M/')

# =======================
# Loading parameters
# =======================
CHUNKSIZE = 10_000_000

# =======================
# PCA Parameters
# =======================

# Number of components for PCA
PCA_N_COMPONENTS = 3

# Number of steps to divide each dimension. len(STEPS_LIST) == len(PCA_N_COMPONENTS)
STEPS_LIST = [64, 32, 16]

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


# =======================
# TMAP Configuration
# =======================

TMAP_NAME = 'representative_cluster'

# =======================
# File Paths
# =======================

INPUT_TMAP_PATH = os.path.join(BASE_DIR, 'data/230M','cluster_representatives.csv')

OUTPUT_TMAP_PATH = os.path.join(BASE_DIR, 'data/230M/', 'maps/')

# =======================
# LSH Configuration
# =======================

PERMUTATIONS = 512 # Number of permutations to be used in MinHashing

TMAP_K = 20 # Number of neighbors 

# =======================
# Layout Configuration
# =======================

TMAP_NODE_SIZE = 1/40

TMAP_POINT_SCALE = 3.0