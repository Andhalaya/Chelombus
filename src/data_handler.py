import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import IncrementalPCA
from src.utils.helper_functions import find_input_type
import numpy as np
import time
import tqdm
 
def get_total_chunks(file_path, chunksize):
    """ Calculate number of chunks based on self.chunksize for tqdm 
    Maybe avoid for files that are too large >150 GB? Takes about ~2 minutes for such size
    You can also just add manually the amount of lines to get an estimation
    """
    print('Preparing tqdm...')
    total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for header
    # total_lines = int(664075400) # Lines for the Enamine_REAL file ~664M compounds
    total_chunks = (int(total_lines) + int(chunksize) - 1) // int(chunksize)
    return total_chunks
    
class DataHandler:
    def __init__(self, file_path, chunksize):
        self.file_path = file_path
        self.chunksize = chunksize
        self.datatype = find_input_type(file_path)

    def get_total_lines(self):
        """Calculate the total number of lines in the file."""
        with open(self.file_path, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)

    def load_data(self):
        total_chunks = get_total_chunks(self.file_path, self.chunksize)

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
        
    def _load_txt_data(self):
        #TODO return data from txt file
        #it should work with load cxsmiles function
        try:
            pass
        except Exception as e:
            raise ValueError(f"Error reading file: {e}")
            
    def _load_cxsmiles_data(self):  # Generator object, to be used in enumerate
        try:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                # Skip the header line
                header = f.readline()  
                
                while True:
                    smiles_and_features = []
                    try:
                        for _ in range(self.chunksize):
                            line = f.readline()
                            if not line:
                                break
                            
                            # Extract the SMILES from the line
                            line_splitted = line.strip().split('\t') # -> List 
                            if line_splitted:  # Ensure it's not an empty line
                                smiles_and_features.append(line_splitted)

                    except Exception as e:
                        raise ValueError(f"Error reading chunk: {e}")

                    if not smiles_and_features:
                        break  # Stop when no more lines

                    yield smiles_and_features  # Yield the chunk of SMILES

        except FileNotFoundError:
            raise ValueError(f"File not found: {self.file_path}")
        except Exception as e:
            raise ValueError(f"Error opening file: {e}")

        
    
    def extract_smiles_and_features(self, data):
            """ 
            Method for extracting smiles and features in pandas. `data` object needs to be pd.DataFrame. 
            Not optimal for very large datasets? So far I've tried with 10M samples and performed well
            input: 
            param: data. pd.Dataframe, .txt or .cxsmiles file 
            output:
            smiles = np.array(batch_size,)
            features = np.array(batch_size, features) 
            """
            try: 
                data.columns = data.columns.str.lower() # for csv file
                try: 
                    smiles_column = data.filter(like='smiles').columns[0]  # Gets the first matching column name
                except: 
                    raise ValueError(f"No SMILES column found in {self.file_path}. Ensure there's a column named 'smiles'.")
                
                smiles_list = data[smiles_column]
                features_list = data.drop(columns=[smiles_column])
            
            except: 
                # for .txt or cxsmiles files
                """
                 This part assumes that the 'smiles' columns in the txt or cxsmile files is in first position
                """
                smiles_list = np.array(data)[:,0] # Return list of all first elements in the list of lists (data) -> list of smiles
                features_list = np.array(data)[:,1:]

                # return smiles_list, features_list
 
            
            return np.array(smiles_list), np.array(features_list)
    
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
    

