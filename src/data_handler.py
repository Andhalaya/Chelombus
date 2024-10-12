from config import DATA_FILE_PATH
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import IncrementalPCA
import numpy as np

class DataHandler:
    def load_data(self, file_path=None, chunksize=None):
        if file_path.endswith('csv'):
                       
            # Load data from a file with optional chunking
            if not file_path:
                file_path = DATA_FILE_PATH
            if not chunksize:
                chunksize = self.chunksize
            try:
                return pd.read_csv(file_path, chunksize=chunksize)
            
            except Exception as e:
                raise ValueError(f"Error reading file: {e}")

        # TODO: Add support for txt  
        # elif file_path.split('.')[-1] == 'txt':
        #     return False 
        
        else: 
            raise ValueError('Unsupported input file. Only .csv files are supported')
    
    
    def extract_smiles_and_features(self, data):
            data.columns = data.columns.str.lower()

            try: 
                smiles_column = data.filter(like='smiles').columns[0]  # Gets the first matching column name
                
            except:
                raise ValueError(f"No SMILES column found in {DATA_FILE_PATH}. Ensure there's a column named 'smiles'.")
                

            smiles_list = data[smiles_column]
            features = data.drop(columns=[smiles_column])

            return np.array(smiles_list), np.array(features)
    
    def get_total_chunks(self, file_path, chunksize):
        total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for header
        total_chunks = (total_lines + chunksize - 1) // chunksize
        return total_chunks
    

    def one_hot_encode(self, features): 
        one_hot_encoder = ColumnTransformer(
            transformers=[
                ('cat', OneHotEncoder(), features)
            ], 
            remainder= 'passthrough' # Keep the rest of columns as is
        )

        oh_features = one_hot_encoder.fit_transform(features)
        return oh_features
    
    def standarize_data(self, data): 
        scaler = StandardScaler()
        data_standardized = scaler.fit_transform(data)

