import hnswlib
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import tmap as tm
from faerun import Faerun
import numpy as np
from sklearn.neighbors import NearestNeighbors
import tmap as tm
from faerun import Faerun
from time import time


file_paths =[]
for i in range(4):
    file_paths.append(f'/home/afloresep/work/chelombus/data/output/batch_data_{i}.csv')


data = pd.concat([pd.read_csv(file) for file in file_paths])

# Extract PCA coordinates (fingerprints)
pca_columns = ['PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']
fingerprints = data[pca_columns].values.astype('float32')  # Convert to float32 for HNSWlib compatibility

# Initialize HNSW index
dim = fingerprints.shape[1]  # Dimension is 5 for PCA vectors
num_elements = len(fingerprints)

# Create the index
p = hnswlib.Index(space='l2', dim=dim)
p.init_index(max_elements=num_elements, ef_construction=200, M=16)

# Add items to the index
p.add_items(fingerprints)

# Perform k-NN query (e.g., k=20)
k = 20
labels, distances = p.knn_query(fingerprints, k=k)

# Construct the edge list
edge_list = []
for i in range(num_elements):
    for j in range(1, k):  # Start from 1 to avoid self-loop
        edge_list.append([i, labels[i][j], distances[i][j]])

# Step 2: Extract PCA feature vectors
pca_columns = ['PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']
fingerprints = data[pca_columns].values

# Step 3: Set up the k-NN model for finding nearest neighbors in 5D space
knn = 20
print('Starting KNN')
knn_search = NearestNeighbors(n_neighbors=knn, metric='manhattan', n_jobs=-1)  # Using 'manhattan' as the metric
knn_search.fit(fingerprints)


print('Creating layout')
# Step 5: Create the Tmap layout configuration
cfg = tm.LayoutConfiguration()
cfg.node_size = 1 / 30
cfg.mmm_repeats = 2
cfg.sl_extra_scaling_steps = 5
cfg.k = 20
cfg.sl_scaling_type = tm.RelativeToAvgLength

print('Generating TMAP layout')
# Step 6: Generate the Tmap layout using the edge list
x_, y_, s, t, gp = tm.layout_from_edge_list(len(fingerprints), edge_list, cfg)

# Plotting
f = Faerun(
    view="front",
    coords=False,
    title="",
    clear_color="#FFFFFF",
)

f.add_scatter(
    "mqn_pca5_TMAP",
    {
        "x": x_,
        "y": y_,
        "c": np.arange(len(fingerprints)),
        "labels": data['smiles'],
    },
    point_scale=3,
    colormap=['rainbow'],
    has_legend=True,
    legend_title=['Standard_value'],
    categorical=[False],
    shader='smoothCircle'
)


f.add_tree("mqn_pca5_TMAP_tree", {"from": s, "to": t}, point_helper="mqn_pca5_TMAP")
f.plot('mqn_pca5_TMAP', template='smiles')
