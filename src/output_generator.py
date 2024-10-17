import pandas as pd
import pickle
import os
import numpy as np
import TDigest

class OutputGenerator():
    def __init__(self):
        # # TODO: Dynamic method for calculating the percentiles for every PCA dimension 
        # digests = {}
        # for i in range(5):
        #     digests[f'{i}_digest'] = TDigest()
    
        self.x_digest = TDigest()
        self.y_digest = TDigest()
        self.z_digest = TDigest()
        
    def _online_batch_update(self, coordinates):
        # Update the digest in batches
        self.x_digest.batch_update(coordinates[:,0])
        self.y_digest.batch_update(coordinates[:,1])
        self.z_digest.batch_update(coordinates[:,2])

    def online_update(self, coordinates):
        # Update the digest sequentially
        for x in range(coordinates):
            self.x_digest.update(coordinates)
        
    def get_percentiles(self, coordinates):
        # Update the digest in batches
        percentiles = []

        self.x_digest.batch_update(coordinates[:,0])
        self.y_digest.batch_update(coordinates[:,1])
        self.z_digest.batch_update(coordinates[:,2])
        
        x_percentile = (self.x_digest.percentile(0.01), self.x_digest.percentile(99.99))
        y_percentile = (self.y_digest.percentile(0.01), self.y_digest.percentile(99.99))
        z_percentile = (self.z_digest.percentile(0.01), self.z_digest.percentile(99.99))

        percentiles.append(x_percentile)
        percentiles.append(y_percentile)
        percentiles.append(z_percentile)

        return percentiles
    
    # Function to map coordinates into [0, 100] steps
    def map_to_grid(coord, percentiles, min_val, step_size=100):

        # Define the percentile ranges for each coordinate (based on 0.01 and 99.99 percentiles)
        x_percentiles = [-23.65769251, 30.68569253]
        y_percentiles = [-16.20908626, 23.04514974]
        z_percentiles = [-10.74670722, 12.98772289]

        # Calculate the step sizes based on the percentile ranges
        x_step_size = (x_percentiles[1] - x_percentiles[0]) / 100
        y_step_size = (y_percentiles[1] - y_percentiles[0]) / 100
        z_step_size = (z_percentiles[1] - z_percentiles[0]) / 100


        # Apply mapping to x, y, and z coordinates
        mapped_x = coord['x'], x_percentiles[0], x_step_size
        mapped_y = coord['y'], y_percentiles[0], y_step_size
        mapped_z = coord['z'], z_percentiles[0], z_step_size

        return np.clip(np.floor((coord - min_val) / step_size), 0, (step_size - 1))

    def save_batch(self, batch_idx, coordinates, smiles_list, features):
        
         batch_data = pd.DataFrame({
             'smiles': smiles_list
         })

         for i in range(len(coordinates[0])):
            batch_data[f'PCA_{i+1}'] = coordinates[:, i]

         batch_data.to_csv(f'data/1M/output/batch_data_{batch_idx}.csv', index=False)