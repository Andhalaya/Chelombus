from config import DATA_FILE_PATH, SMILES_COLUMN_NAME

class DataHanlder():

    def __init__(self):
        self.data = None
        self.smiles_column = SMILES_COLUMN_NAME
        

    def load_data(self, file_path=DATA_FILE_PATH): # read data and detect file format automatically
        pass

    def extract_smiles(): # accept chunk as an argument
        pass 
        
    def extract_features(self): #dynamically extract all columns excep the SMILES 
        # Should accept chunk as an argument
        self.features = self.data.drop(columns=[self.simles_column])
        pass
        # if no additional features are present, pipeline should function, fetures variable can be an empty DF or None