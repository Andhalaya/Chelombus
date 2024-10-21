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


    def save_batch(self, batch_idx, coordinates, smiles_list, features):
         """
         Save the batch output as CSV file
         """
         batch_data = pd.DataFrame({'smiles': smiles_list})

         for i in range(len(coordinates[0])):
            batch_data[f'PCA_{i+1}'] = coordinates[:, i]
        
            output_path = os.path.join(OUTPUT_FILE_PATH, f'batch_data_{batch_idx}.csv')  

            os.makedirs(OUTPUT_FILE_PATH, exist_ok=True)

            batch_data.to_csv(output_path, index=False)


    def _round_to_step(self,coordinate:float, min_value:float, max_value:float, step_size:float):
           
            # Map the coordinate to its closest value in the steps
            
            if coordinate < min_value: 
                    return min_value
            elif coordinate > max_value:
                    return max_value
            else: # min_value + Number of steps x Step Size
                    return (min_value) + (step_size)*round((coordinate - min_value)/step_size)
                    

    def fit_coordinates(self, output: str, percentiles:list, grid:list):
            #TODO: Handle more dimensions than 3 

            """
            Fit the PCA coordinates of every point x=[PCA_1, PCA_2, PCA_3] to the grid of the cub. 
            output: previously generated csv with smiles,PCA_1,PCA_2, PCA_3 columns
            percentiles: percentiles to be used as range for every dimension
            dimensions: number of steps in every range
             
            e.g. percentiles [(-10, 10), (0,10), (-20, 10)] and steps = [1, 1, 5]
            will create a cube where
            values in the x-axis will take coordinates from -10,10 in steps of 10: -10, -9,..., 9, 10
            values in the y-axis will take coordinates from 0,10 in steps of 1:  0,1, 2, ... , 9, 10 
            values in the z-axis will take coordinates from -20, 10 in steps of: -20, -15,..., 15, 20
            """

            step_PCA_1 = (percentiles[0][1] - percentiles[0][0])/ grid[0] 
            step_PCA_2 = (percentiles[1][1] - percentiles[1][0])/ grid[1] 
            step_PCA_3 = (percentiles[2][1] - percentiles[2][0])/ grid[2] 
            
            df_output = pd.read_csv(os.path.join(OUTPUT_FILE_PATH, output))
            
            df_output['PCA_1'] = df_output['PCA_1'].apply(lambda x: self._round_to_step(x, percentiles[0][0], percentiles[0][1], step_PCA_1))
            df_output['PCA_2'] = df_output['PCA_2'].apply(lambda x: self._round_to_step(x, percentiles[1][0], percentiles[1][1], step_PCA_2))
            df_output['PCA_3'] = df_output['PCA_3'].apply(lambda x: self._round_to_step(x, percentiles[2][0], percentiles[2][1], step_PCA_3))

            output_path = os.path.join(OUTPUT_FILE_PATH, output)
            df_output.to_csv(output_path, index=False)

    def fit_coord_multidimensional(self, output:str , percentiles: list):
            """
            Generalize fit_coordinates to handle more than 3 dimensions.
            This function will handle any number of dimensions.

            :param output: The CSV file containing the PCA coordinates.
            :param percentiles: A list of percentile ranges for each dimension.
            :param grid: A list of steps for each dimension.
            """
            if len(percentiles) != len(self.steps):
                 raise ValueError("Percentiles and grid must have same number of dimensions")
            
            df_output = pd.read_csv(os.path.join(OUTPUT_FILE_PATH, output))

            for dim in range(len(percentiles)):
                step = (percentiles[dim][1] - percentiles[dim][0]) / self.steps[dim]
                pca_column = f'PCA_{dim + 1}'

                # Apply rounding for each dimension
            df_output[pca_column] = df_output[pca_column].apply(
            lambda x: self._round_to_step(x, percentiles[dim][0], percentiles[dim][1], step)
            )
