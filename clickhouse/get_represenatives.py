import time
import pandas as pd
import clickhouse_connect
import numpy as np
from scipy.spatial.distance import cdist
from dynamic_tmap.src.utils.helper_functions import ProgressTracker, format_time


# Connect to DB 
client = clickhouse_connect.get_client(host='localhost', port=8123)

# Execute the query and load data into a DataFrame
query = "SELECT DISTINCT cluster_id FROM clustered_enamine"
query_all_clusters = client.query(query)

all_clusters = query_all_clusters.result_rows

# Initialize empty lists to store representative molecule data
representative_medians = []

s = time.time()
cluster_number = 0
# Iterate through the clusters
with ProgressTracker(description="Clustering representatives", total_steps=len(all_clusters)) as tracker:
    for cluster_id in all_clusters:
        s = time.time()
        cluster_number += 1 
        print(f"Loading cluster: {cluster_id[0]}. {cluster_number}/{len(all_clusters)}")
        # Load the dataframe for the current cluster
        df = pd.DataFrame(
            client.query(f"SELECT * FROM clustered_enamine WHERE cluster_id == {cluster_id[0]}").result_rows,
            columns=["smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id"]
        )
        
        # Median-based approach
        median_coords = df[['PCA_1', 'PCA_2', 'PCA_3']].median().to_numpy()
        median_index = np.argmin(  # Find the molecule closest to the median in PCA space
            cdist([median_coords], df[['PCA_1', 'PCA_2', 'PCA_3']].to_numpy(), metric='euclidean')
        )
        median_molecule = df.iloc[median_index]  # Retrieve the median-based representative molecule
        
        # Save the median-based representative molecule details
        representative_medians.append({
            "smiles": median_molecule["smiles"],
            "PCA_1": median_molecule["PCA_1"],
            "PCA_2": median_molecule["PCA_2"],
            "PCA_3": median_molecule["PCA_3"],
            "cluster_id": median_molecule["cluster_id"]
        })
        e = time.time()

        tracker.update_progress()
        # print(f"Cluster {cluster_id} done: Elapsed time {format_time(e - s)}")

# Convert the lists of representative molecules into new dataframes
representatives_median_df = pd.DataFrame(representative_medians)

representatives_median_df.to_csv('/mnt/10tb_hdd/clustered_enamine_database/enamine_representatives.csv', index=False)

e = time.time()
print(f"Elapsed time: {format_time(e - s)} hours")