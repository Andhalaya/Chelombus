import pandas as pd
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import umap
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.express as px


from rdkit import Chem
from rdkit.Chem import rdMolDescriptors                                                           
import numpy as np

import os
import time
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Step 1: Load CSV Data and Sample a Subset
# Replace the filenames with your actual CSV files
data = pd.concat([pd.read_csv(file) for file in ["/home/afloresep/work/chelombus/data/output/batch_data_0.csv", "/home/afloresep/work/chelombus/data/output/batch_data_1.csv", '/home/afloresep/work/chelombus/data/output/batch_data_2.csv', '/home/afloresep/work/chelombus/data/output/batch_data_3.csv']])  # Adjust for your filenames



def smiles_to_fingerprint(smiles): 
    fingerprint = rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smiles))
    return np.array(fingerprint)

# Apply fingerprint generation to the sampled data
data['fingerprint'] = data['smiles'].apply(smiles_to_fingerprint)

# Step 3: Dimensionality Reduction using UMAP
fingerprint_matrix = np.vstack(data['fingerprint'].values)

# Use UMAP to reduce dimensionality to 2D or 3D
umap_reducer = umap.UMAP(n_components=2)  # Change to n_components=3 for 3D
umap_embedding = umap_reducer.fit_transform(fingerprint_matrix)

# Step 4: Visualization
# 2D Scatter Plot
plt.figure(figsize=(10, 7))
plt.scatter(data['x'], data['y'], c=umap_embedding[:, 0], cmap='viridis', alpha=0.6)
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('2D Visualization of SMILES Similarity')
plt.colorbar(label='UMAP Component')
plt.savefig('2D-similarity.png')

# 3D Scatter Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(data['x'], data['y'], data['z'], c=umap_embedding[:, 0], cmap='viridis', alpha=0.6)
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')
ax.set_title('3D Visualization of SMILES Similarity')
plt.colorbar(sc, label='UMAP Component')
plt.savefig('3D visualization.png')


# Step 5: Interactive Plot using Plotly
fig = px.scatter_3d(
    data,
    x='x',
    y='y',
    z='z',
    color=umap_embedding[:, 0],
    hover_data=['smiles'],
    title='Interactive 3D Visualization of SMILES Similarity'
)
fig.update_traces(marker=dict(size=5, opacity=0.7))
fig.show()

# Note: Ensure you have all required libraries installed
# Install RDKit: conda install -c conda-forge rdkit
# Install UMAP: pip install umap-learn
# Install Plotly: pip install plotly
# Install Matplotlib, Pandas: pip install matplotlib pandas

# Adjust the filenames and sampling parameters as per your data.