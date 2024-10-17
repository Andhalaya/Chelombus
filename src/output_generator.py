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
             'smiles': smiles_list
         })

         for i in range(len(coordinates[0])):
            batch_data[f'PCA_{i+1}'] = coordinates[:, i]

         batch_data.to_csv(f'data/output/batch_data_{batch_idx}.csv', index=False)