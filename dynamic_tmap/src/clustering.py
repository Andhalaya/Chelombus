import dask.dataframe as dd
import math
import os

# Read the parquet files into a Dask DataFrame
ddf = dd.read_parquet("/mnt/samsung_2tb/mixed_data/output/pca_transformed_results/output/*.parquet")

# Sort the entire dataset by PCA_1
# Using set_index with sorted=True ensures that the data is globally sorted by PCA_1
ddf = ddf.set_index("PCA_1")

# Reset the index to get a global integer index from 0 to (total_rows - 1)
# PCA_1 now becomes a column again, and the new index is a simple RangeIndex across all rows
ddf = ddf.reset_index(drop=False)

# Define rows per file and calculate the number of files
rows_per_file = 1_340_000
total_rows = ddf.shape[0].compute()
num_files = math.ceil(total_rows / rows_per_file)

# Create the output directory
output_dir = "/mnt/samsung_2tb/dask_sorted_chunks/"
os.makedirs(output_dir, exist_ok=True)

# Iterate over the desired number of files, slicing by row index range
for i in range(num_files):
    start_idx = i * rows_per_file
    end_idx = (i + 1) * rows_per_file

    # Filter based on the global integer index
    subset = ddf[(ddf.index >= start_idx) & (ddf.index < end_idx)]
    
    # Since we're filtering a Dask DataFrame by row-based conditions, we need to compute before writing
    subset_df = subset.compute()

    # Write to Parquet file
    output_filename = os.path.join(output_dir, f"sorted_subset_{i+1:04d}.parquet")
    subset_df.to_parquet(output_filename)
    
    print(f"Saved file {i+1}/{num_files}: {output_filename}")
