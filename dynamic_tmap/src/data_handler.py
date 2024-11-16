import gc 
import logging
import os 
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import IncrementalPCA
from src.utils.helper_functions import find_input_type
from memory_profiler import profile
import numpy as np
import time
import tqdm
from src.fingerprint_calculator import FingerprintCalculator
    
class DataHandler:

    def __init__(self, file_path, chunksize):
        self.file_path = file_path
        self.chunksize = chunksize
        self.datatype = find_input_type(file_path)

    @staticmethod
    def get_total_chunks(file_path, chunksize):
        """
        Calculate number of chunks based on self.chunksize for tqdm 
        Maybe avoid for files that are too large >150 GB? Takes about ~2 minutes for such size
        You can also just add manually the amount of lines to get an estimation
        """
        print('Preparing tqdm...')
        total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for header
        # total_lines = int(664075400) # Lines for the Enamine_REAL file ~664M compounds
        total_chunks = (int(total_lines) + int(chunksize) - 1) // int(chunksize)
        return total_chunks

    def get_total_lines(self):
        """Calculate the total number of lines in the file."""
        with open(self.file_path, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)

    def load_data(self):
        """Loads data in chunks depending on filetype. 
        Output:
          - Generator function according to datatype
          - Number of chunks (for tqdm)
        """
        total_chunks = self.get_total_chunks(self.file_path, self.chunksize)

        # Dynamically dispatch the right method based on datatype
        if self.datatype == 'csv':
            return self._load_csv_data(), total_chunks
        
        elif self.datatype == 'cxsmiles' :
            return self._load_cxsmiles_data(), total_chunks
        else:
            raise ValueError(f"Unsupported file type: {self.datatype}")

    def _load_csv_data(self): # Generator object
        try:
            return pd.read_csv(self.file_path, chunksize=self.chunksize)

        except Exception as e:
            raise ValueError(f"Error reading file: {e}")

    def _load_cxsmiles_data(self):
        smiles_col_index = 1  # Index for the 'smiles' column (0-based), first column was hash value after sorting
        try:
            for chunk in pd.read_csv(
                self.file_path,
                sep='\t',
                header=None,         # No header row
                usecols=[smiles_col_index],
                chunksize=self.chunksize,
                dtype=str,
                engine='c',
                encoding='utf-8'
            ):
                # smiles_list = chunk[smiles_col].tolist() -> for when data has column names
                smiles_list = chunk.iloc[:, 0].tolist()  # Access the first column of the chunk
                yield smiles_list
        except Exception as e:
            raise ValueError(f"Error loading data: {e}")


    def process_chunk(self, idx, chunk, output_dir):
        """
        Process a single chunk of data by calculating fingerprints and saving them to a parquet file
        """
        try:
            # Check if chunk already exists
            fp_chunk_path = os.path.join(output_dir, f'batch_parquet/fingerprints_chunk_{idx}.parquet')
            if os.path.exists(fp_chunk_path):
                logging.info(f'Chunk {idx} already processed, skipping.')
                return            

            # Extract smiles and features from chunk
            fp_calculator = FingerprintCalculator(chunk, 'mqn')

            # Calculate fingerprints
            fingerprints = fp_calculator.calculate_fingerprints()
    
            # Ensure output directories exist
            os.makedirs(os.path.join(output_dir, 'batch_parquet'), exist_ok=True)
            os.makedirs(os.path.join(output_dir, 'output'), exist_ok=True)

            # Create dataframe with smiles list
            smiles_dataframe= pd.DataFrame({
                'smiles': smiles_list, 
            })

            del smiles_list

            # Create dataframe with fingerprints values 
            fingerprint_df = pd.DataFrame(fingerprints.tolist(), columns = [f'fp_{i+1}' for i in range(42)])
            del fingerprints

            # Concat both df. 
            chunk_dataframe = pd.concat([smiles_dataframe, fingerprint_df], axis=1)
            del smiles_dataframe, fingerprint_df 

            # Save to parquet dataframe
            chunk_dataframe.to_parquet(fp_chunk_path, index=False)

            del chunk_dataframe
            gc.collect() # Collect garbage. Not sure if it makes a difference but just in case

        except Exception as e:
            logging.error(f"Error processing chunk {idx}: {e}", exc_info=True)