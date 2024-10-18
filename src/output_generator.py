import pandas as pd
import pickle
import os
from config import OUTPUT_FILE_PATH, CHUNKSIZE

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
        
            output_path = os.path.join(OUTPUT_FILE_PATH, f'batch_data_{batch_idx}.csv')  

            batch_data.to_csv(output_path, index=False)


    def _round_to_step(self,coordinate:float, min_value:float, max_value:float, step_size:float):
           
            # Map the coordinate to its closest value in the steps
            
            if coordinate < min_value: 
                    return min_value
            elif coordinate > max_value:
                    return max_value
            else: # min_value + Number of steps x Step Size
                    mapped_coordinate = (min_value) + (step_size)*round((coordinate - min_value)/step_size)
                    return mapped_coordinate

    def fit_coordinates(self, output: str, percentiles:list):

            step_PCA_1 = (percentiles[0][1] - percentiles[0][0])/ 10 
            step_PCA_2 = (percentiles[1][1] - percentiles[1][0])/ 10
            step_PCA_3 = (percentiles[2][1] - percentiles[2][0])/ 10
            
            df_output = pd.read_csv(os.path.join(OUTPUT_FILE_PATH, output))
            
            df_output['PCA_1'] = df_output['PCA_1'].apply(lambda x: self._round_to_step(x, percentiles[0][0], percentiles[0][1], step_PCA_1))
            df_output['PCA_2'] = df_output['PCA_2'].apply(lambda x: self._round_to_step(x, percentiles[1][0], percentiles[1][1], step_PCA_2))
            df_output['PCA_3'] = df_output['PCA_3'].apply(lambda x: self._round_to_step(x, percentiles[2][0], percentiles[2][1], step_PCA_3))

            df_output.to_csv(os.path.join(OUTPUT_FILE_PATH, output), index=False)
