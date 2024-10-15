import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import IncrementalPCA
from src.utils.helper_functions import find_input_type
import numpy as np
import time

def get_total_chunks(file_path, chunksize):
    """ Calculate number of chunks based on self.chunksize. For tqdm 
    avoid for files that are too large >150 GB? Takes about ~2 minutes for such size
    """
    total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for header
    total_chunks = (total_lines + chunksize - 1) // chunksize
    return total_chunks

class DataHandler:
    def __init__(self, file_path, chunksize):
        self.file_path = file_path
        self.chunksize = chunksize

    def load_data(self):      
        if find_input_type(self.file_path) == 'csv':       
            try:
                total_chunks = get_total_chunks(self.file_path, self.chunksize)
                
                return pd.read_csv(self.file_path, chunksize=self.chunksize), total_chunks
            
            except Exception as e:
                raise ValueError(f"Error reading file: {e}")

        elif find_input_type(self.file_path) == '.txt':
            # TODO: Add support for txt  
            pass 
        
        elif self.file_path.split('.')[-1] == '.cxsmiles':
            try:
                with open(self.file_path,'r') as f:
                    pass
            except Exception as e:
                raise ValueError(f"Error reading file: {e}")

        else: 
            raise ValueError('Unsupported input file. Only .csv, .txt. and .CXSSmiles files are supported')
    
    
    def extract_smiles_and_features(self, data):
            """ 
            Method for extracting smiles and features in pandas. `data` object needs to be pd.DataFrame. 
            Not optimal for very large datasets? So far I've tried with 10M samples and performed well
            """
            data.columns = data.columns.str.lower()

            try: 
                smiles_column = data.filter(like='smiles').columns[0]  # Gets the first matching column name
                
            except:
                raise ValueError(f"No SMILES column found in {self.file_path}. Ensure there's a column named 'smiles'.")
                

            smiles_list = data[smiles_column]
            features = data.drop(columns=[smiles_column])

            return np.array(smiles_list), np.array(features)
    



    def one_hot_encode(self, features): 
        """ 
        TODO Idea is to one hot encode categorical features e.g. 'target value'
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
    

