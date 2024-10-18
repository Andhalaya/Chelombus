import pandas as pd
import pickle
import os

def get_percentiles(x_digest, y_digest, z_digest):
  
    # Update the digest in batches
    percentiles = []
  
    x_percentile = (x_digest.percentile(0.01), x_digest.percentile(99.99))
    y_percentile = (y_digest.percentile(0.01), y_digest.percentile(99.99))
    z_percentile = (z_digest.percentile(0.01), z_digest.percentile(99.99))
    percentiles.append(x_percentile)
    percentiles.append(y_percentile)
    percentiles.append(z_percentile)

    return percentiles 


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

         batch_data.to_csv(f'data/10M/output/batch_data_{batch_idx}.csv', index=False)

         