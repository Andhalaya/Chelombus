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
        
            output_path = os.path.join(OUTPUT_FILE_PATH, f'output/batch_data_{batch_idx}.csv')  

            os.makedirs(os.path.join(OUTPUT_FILE_PATH, 'output'), exist_ok=True)

            batch_data.to_csv(output_path, index=True)


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
            Generalize fit_coordinates to handle more than 3 dimensions.
            This function will handle any number of dimensions.

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

            df_output.to_csv(output_path, index=False)

