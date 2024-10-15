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
from sklearn.cluster import KMeans
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors        
from sklearn.cluster import DBSCAN
import os
import pandas as pd
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
import umap
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.express as px
from sklearn.cluster import KMeans
import os

# Step 1: Load CSV Data and Sample a Subset
# Replace the filenames with your actual CSV files
data = pd.concat([pd.read_csv(file) for file in ["/home/afloresep/work/chelombus/data/output/batch_data_0.csv",
                                                  "/home/afloresep/work/chelombus/data/output/batch_data_1.csv"]])  # Adjust for your filenames

# Randomly sample 5000 points for visualization
sampled_data = data.sample(5000)

# Step 2: Generate Molecular Fingerprints

def smiles_to_fingerprint(smiles): 
    fingerprint = rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smiles))
    return np.array(fingerprint)


# Apply fingerprint generation to the sampled data
sampled_data['fingerprint'] = sampled_data['smiles'].apply(smiles_to_fingerprint)

# Step 3: Dimensionality Reduction using UMAP
fingerprint_matrix = np.vstack(sampled_data['fingerprint'].values)

# Use UMAP to reduce dimensionality to 2D or 3D
umap_reducer = umap.UMAP(n_components=3)  # Change to n_components=3 for 3D
umap_embedding = umap_reducer.fit_transform(fingerprint_matrix)

# Add UMAP 3D coordinates to sampled_data
sampled_data[['umap_x', 'umap_y', 'umap_z']] = umap_embedding


# Step 4: Clustering Based on 3D Coordinates
# Use KMeans to cluster points, targeting approximately 10 points per cluster
num_clusters = len(sampled_data) // 10
kmeans = KMeans(n_clusters=num_clusters, random_state=42).fit(umap_embedding)
sampled_data['cluster'] = kmeans.labels_

# Save clustering information to a text file
with open('clusters.txt', 'w') as f:
    for cluster_id in sorted(sampled_data['cluster'].unique()):
        f.write(f'Cluster {cluster_id}')
        cluster_points = sampled_data[sampled_data['cluster'] == cluster_id]
        for _, row in cluster_points.iterrows():
            f.write(f"SMILES: {row['smiles']}, Coordinates: ({row['umap_x']}, {row['umap_y']}, {row['umap_z']})\n")
        f.write('\n')

# Step 5: Visualization
# 2D Scatter Plot
plt.figure(figsize=(10, 7))
plt.scatter(sampled_data['x'], sampled_data['y'], c=umap_embedding[:, 0], cmap='viridis', alpha=0.6)
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('2D Visualization of SMILES Similarity')
plt.colorbar(label='UMAP Component')
# plt.savefig()

# 3D Scatter Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(sampled_data['x'], sampled_data['y'], sampled_data['z'], c=umap_embedding[:, 0], cmap='viridis', alpha=0.6)
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')
ax.set_title('3D Visualization of SMILES Similarity')
plt.colorbar(sc, label='UMAP Component')
# plt.savefig()

# Step 6: Interactive Plot using Plotly
fig = px.scatter_3d(
    sampled_data,
    x='x',
    y='y',
    z='z',
    color=umap_embedding[:, 0],
    hover_data=['smiles'],
    title='Interactive 3D Visualization of SMILES Similarity'
)
fig.update_traces(marker=dict(size=5, opacity=0.7))
# fig.show()


# Step 7: Visualize Each Cluster
output_dir = "cluster_visualizations"
os.makedirs(output_dir, exist_ok=True)

for cluster_id in sorted(sampled_data['cluster'].unique()):
    cluster_points = sampled_data[sampled_data['cluster'] == cluster_id]
    smiles_list = cluster_points['smiles'].tolist()
    mols = [Chem.MolFromSmiles(smile) for smile in smiles_list]
    
    # Generate 2D coordinates for each molecule
    for mol in mols:
        if mol is not None:
            AllChem.Compute2DCoords(mol)
    
    # Create a subplot for each cluster
    n_mols = len(mols)
    n_cols = min(6, n_mols)
    n_rows = (n_mols + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3))
    axes = axes.flatten() if n_mols > 1 else [axes]
    fig.suptitle(f"Cluster {cluster_id} - 2D Structures of SMILES", fontsize=16)
    
    for i, (mol, ax) in enumerate(zip(mols, axes)):
        if mol is not None:
            img = Draw.MolToImage(mol)
            ax.imshow(img)
            ax.axis('off')
            ax.set_title(f"Structure {i+1}", fontsize=10)
    
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(output_dir, f"cluster_{cluster_id}.png"))
    plt.close(fig)