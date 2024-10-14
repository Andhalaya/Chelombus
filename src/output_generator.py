import pandas as pd
import pickle
import os

class OutputGenerator():
    def __init__(self):
        pass 

    def compile_output(self, coordinates, smiles_list, features): 
        pass

    def save_batch(self, batch_idx, coordinates, smiles_list, features):
         batch_data = pd.DataFrame({
             'smiles': smiles_list,  
             'x': coordinates[:, 0], 
             'y': coordinates[:, 1], 
             'z': coordinates[:, 2], 
         })
         batch_data.to_csv(f'data/output/batch_data_{batch_idx}.csv', index=False)