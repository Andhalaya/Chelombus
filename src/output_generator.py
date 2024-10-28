import numpy as np
import pandas as pd
import pickle
import os
from config import OUTPUT_FILE_PATH, CHUNKSIZE

def get_percentiles(x_digest, y_digest, z_digest):

    """
    Get percentiles for x, y, z dimensions.
    """
    percentiles = [(x_digest.percentile(0.01), x_digest.percentile(99.99)),
                   (y_digest.percentile(0.01), y_digest.percentile(99.99)),
                   (z_digest.percentile(0.01), z_digest.percentile(99.99))]
    return percentiles 

class OutputGenerator():

    @property
    def steps(self):
        """Getter for steps."""
        if self._steps is None:
            raise ValueError("Steps have not been set yet")
        return self._steps
    

    @steps.setter
    def steps(self, steps_list):
        """Setter for steps with validation."""
        if any(step <= 0 for step in steps_list):
            raise ValueError("All steps must be positive values.")
        self._steps = steps_list

    def compile_output(self, coordinates, smiles_list, features):
        """
        Placeholder for compiling output. Add specific implementation here.
        """
        pass


    def save_batch_csv(self, batch_idx, coordinates, smiles_list, features, output_dir):
         """
         Save the batch output as CSV file
         This will generate a csv file for every batch. 
         Not recommended for very large files as will increase I/O
         """
         batch_data = pd.DataFrame({'smiles': smiles_list})

         for i in range(len(coordinates[0])):
            batch_data[f'PCA_{i+1}'] = coordinates[:, i]
        
            output_path = os.path.join(output_dir, f'output/batch_data_{batch_idx}.csv')  

            os.makedirs(os.path.join(output_dir, 'output'), exist_ok=True)

            batch_data.to_csv(output_path, index=True)

    def save_batch_parquet(self, coordinates: np.array, smiles_list:list, features:list, output_dir:str): 
        """
        Similar to save_batch_csv
        This will instead save every batch -coordinates, smiles...- into a single parquet file 
        by merging them every time a chunk is loaded. 
        """
        parquet_path= os.path.join(output_dir, 'output_dataframe.parquet')
        os.makedirs(os.path.join(output_dir, 'output'), exist_ok= True)

        # New dataframe for the chunk that is going to be concatenated to a single file 
        batch_data = pd.DataFrame({'smiles': smiles_list})

        # Append PCA coordinates
        for i in range(len(coordinates[0])):
           batch_data[f'PCA_{i+1}'] = coordinates[:, i]

        # TODO: Add features to DataFrame using python dict instead of list. 
        # for feature_name, feature_values in features.items():
            #  batch_data[feature_name] = feature_values 

        if os.path.exists(parquet_path):
             # Load the existing Parquet file 
             existing_parquet= pd.read_parquet(parquet_path)
             # Append the new batch to the existing parquet file
             concatenated_parquet= pd.concat([existing_parquet, batch_data], ignore_index=True) 
        else: 
             # If parquet file doesn't exit, batch is new 
             concatenated_parquet = batch_data

        concatenated_parquet.to_parquet(parquet_path, engine="pyarrow")

    def _round_to_step(self,coordinate:float, min_value:float, max_value:float, step_size:float):
           
            # Map the coordinate to its closest value in the steps
            
            if coordinate < min_value: 
                    return min_value
            elif coordinate > max_value:
                    return max_value
            else: # min_value + Number of steps x Step Size
                    return (min_value) + (step_size)*round((coordinate - min_value)/step_size)
                  

    def fit_coord_multidimensional(self, output:str , percentiles: list):
            """
            Generalize fit_coordinates that handles more than 3 dimensions.

            :param output: The CSV file containing the PCA coordinates.
            :param percentiles: A list of percentile ranges for each dimension.
            :param grid: A list of steps for each dimension.
            """
            if len(percentiles) != len(self.steps):
                 raise ValueError("Percentiles and grid must have same number of dimensions")
            
            df_output = pd.read_csv(os.path.join(OUTPUT_FILE_PATH, ('output/'+output)))

            for dim in range(len(percentiles)):
                step = (percentiles[dim][1] - percentiles[dim][0]) / self.steps[dim]
                pca_column = f'PCA_{dim + 1}'
                # Apply rounding for each dimension
                df_output[pca_column] = df_output[pca_column].apply(
                lambda x: self._round_to_step(x, percentiles[dim][0], percentiles[dim][1], step)
                )
            
            output_path = os.path.join(OUTPUT_FILE_PATH, ('output/'+output))

            df_output.to_csv(output_path, index=False)# Assuming percentiles is a dictionary with keys 'PCA_1', 'PCA_2', ..., mapping to sorted lists for each PCA column