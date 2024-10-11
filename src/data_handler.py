from config import DATA_FILE_PATH, CHUNKSIZE
import pandas as pd

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

            return smiles_list, features
    
