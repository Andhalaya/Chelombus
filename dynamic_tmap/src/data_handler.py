import gc 
import logging
import os 
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from src.utils.helper_functions import find_input_type
from src.fingerprint_calculator import FingerprintCalculator

def get_total_chunks(file_path, chunksize):
    """ Calculate number of chunks based on self.chunksize for tqdm 
    Maybe avoid for files that are too large >150 GB? Takes about ~2 minutes for such size
    You can also just add manually the amount of lines to get an estimation
    """
    print('Preparing tqdm...')
    total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for header
    total_chunks = (int(total_lines) + int(chunksize) - 1) // int(chunksize)
    return total_chunks
    
class DataHandler:
    def __init__(self, 
                 file_path: str, 
                 chunksize: str, 
                 smiles_col_index: int = 0):
        self.file_path = file_path
        self.chunksize = chunksize
        self.datatype = find_input_type(file_path)
        self.smiles_col_index = smiles_col_index

    def get_total_lines(self):
        """Calculate the total number of lines in the file."""
        with open(self.file_path, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)

    def load_data(self):
        total_chunks = get_total_chunks(self.file_path, self.chunksize)

        # Dynamically dispatch the right method based on datatype
        if self.datatype == 'csv':
            return self._load_csv_data(), total_chunks
        
        elif self.datatype == 'cxsmiles' or self.datatype == 'txt' :
            return self._load_tabular_data(), total_chunks
        else:
            raise ValueError(f"Unsupported file type: {self.datatype}")

    def _load_csv_data(self): # Generator object
        try:
            return pd.read_csv(self.file_path, chunksize=self.chunksize)

        except Exception as e:
            raise ValueError(f"Error reading file: {e}")
        
    def _load_tabular_data(self):
        smiles_col_index = self.smiles_col_index  # Index for the 'smiles' column (0-based)
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

    def one_hot_encode(self, features): 
        """ 
        #TODO Idea is to one hot encode categorical features e.g. 'target value'
        and use it in PCA 
        """
        one_hot_encoder = ColumnTransformer(
            transformers=[
                ('cat', OneHotEncoder(), features)
            ], 
            remainder= 'passthrough' # Keep the rest of columns as is
        )

        oh_features = one_hot_encoder.fit_transform(features)
        return oh_features
    
    # @profile
    def process_chunk(self, 
                      idx, 
                      chunk, 
                      output_dir, 
                      fingerprint: str):
        """
        Process a single chunk of data by calculating fingerprints and saving them to a parquet file
        """
        # TODO: Remove idx and checking if chunk already exists. The reason is there's conflict when dealing with loading from differnt 
        # files. After done with file 1 it will name again chunk_0 from file 2 but will be skipped because chunk_0 from file_1 already exists

        try:
            # Check if chunk already exists
            fp_chunk_path = os.path.join(output_dir, f'batch_parquet/fingerprints_chunk_{idx}.parquet')
            if os.path.exists(fp_chunk_path):
                logging.info(f'Chunk {idx} already processed, skipping.')
                return            

            # Extract smiles and features from chunk
            smiles_list = chunk 

            fp_calculator = FingerprintCalculator(smiles_list, fingerprint_type=fingerprint)

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

            #TODO: Return dataframe instead of saving? Probably the best option
            # Save to parquet dataframe
            chunk_dataframe.to_parquet(fp_chunk_path, index=False)

            del chunk_dataframe
            gc.collect() # Collect garbage. Not sure if it makes a difference but just in case

        except Exception as e:
            logging.error(f"Error processing chunk {idx}: {e}", exc_info=True)

