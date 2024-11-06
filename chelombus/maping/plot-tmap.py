import html
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
from time import time
from tqdm import tqdm
import tmap as tm
import pandas as pd
import numpy as np
from timeit import default_timer as timer
from faerun import Faerun
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from mapchiral.mapchiral import encode as mapc_enc
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import pickle
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import html


def calculate_threshold(data):
    # Function to calculate threshold using IQR method
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    threshold = q3 + 1.5 * iqr
    return threshold

def calculate_molecular_properties(smiles):
    """
    Calculate molecular properties using RDKit.
    
    Args:
    smiles (str): SMILES string of the molecule.
    
    Returns:
    tuple: (HAC, fraction_aromatic_atoms, number_of_rings, clogP, fraction_Csp3)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    hac = mol.GetNumHeavyAtoms()
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    fraction_aromatic_atoms = num_aromatic_atoms / hac if hac > 0 else 0
    number_of_rings = rdMolDescriptors.CalcNumRings(mol)
    clogP = Descriptors.MolLogP(mol)
    fraction_Csp3 = Descriptors.FractionCSP3(mol)
    
    return (hac, fraction_aromatic_atoms, number_of_rings, clogP, fraction_Csp3)


def plot_faerun(x, y, s, t, df):


    """
    Plot the data using Faerun.
    
    Args:
    x (list): X coordinates.
    y (list): Y coordinates.
    s (list): Source nodes for tree plot.
    t (list): Target nodes for tree plot.
    df (pandas.DataFrame): DataFrame with target data.
    """
    f = Faerun(view="front", 
               coords=False, 
               title= "Representatives Cluster", 
               clear_color="#FFFFFF")

    labels = []
    hac_data = []
    frac_aromatic_data = []
    num_rings_data = []
    clogp_data = []
    frac_csp3_data = []

    for i, row in df.iterrows():
        target_name = str(row["cluster_id"]).strip()  # Convert to string and remove leading/trailing whitespace
        target_name = html.escape(target_name)  # Escape special characters
        
        if not target_name:
            target_name = "N/A"  # Provide a default value if empty
        
        labels.append(
            row['smiles']
            + '__'
            + f'<a target="_blank" href="">{target_name}</a><br>'
        )
        
        # Calculate molecular properties
        properties = calculate_molecular_properties(row['smiles'])
        if properties:
            hac, frac_aromatic, num_rings, clogp, frac_csp3 = properties
            hac_data.append(hac)
            frac_aromatic_data.append(frac_aromatic)
            num_rings_data.append(num_rings)
            clogp_data.append(clogp)
            frac_csp3_data.append(frac_csp3)
        else:
            # Handle invalid SMILES
            hac_data.append(None)
            frac_aromatic_data.append(None)
            num_rings_data.append(None)
            clogp_data.append(None)
            frac_csp3_data.append(None)
   
    # Calculate threshold for hac_data using IQR
    hac_threshold = calculate_threshold(hac_data)
    frac_threshold = calculate_threshold(frac_aromatic_data)
    rings_threshold = calculate_threshold(num_rings_data)
    clogp_threshold = calculate_threshold(clogp_data)
    csp3_threshold = calculate_threshold(frac_csp3_data)

    # Function to apply thresholds and return filtered data as separate lists
    def apply_thresholds(hac_data, frac_aromatic_data, num_rings_data, clogp_data, frac_csp3_data):
        filtered_hac = []
        filtered_frac_aromatic = []
        filtered_num_rings = []
        filtered_clogp = []
        filtered_frac_csp3 = []

        # Iterate through all data points and apply thresholds
        for hac, frac, rings, clogp, csp3 in zip(hac_data, frac_aromatic_data, num_rings_data, clogp_data, frac_csp3_data):
            if hac <= hac_threshold and frac <= frac_threshold and rings <= rings_threshold and clogp <= clogp_threshold and csp3 <= csp3_threshold:
                filtered_hac.append(hac)
                filtered_frac_aromatic.append(frac)
                filtered_num_rings.append(rings)
                filtered_clogp.append(clogp)
                filtered_frac_csp3.append(csp3)

        return filtered_hac, filtered_frac_aromatic, filtered_num_rings, filtered_clogp, filtered_frac_csp3

    filtered_hac,filtered_frac_aromatic, filtered_num_rings, filtered_clogp, filtered_frac_csp3 = apply_thresholds(hac_data, frac_aromatic_data, num_rings_data, clogp_data, frac_csp3_data)

    # Add scatter plot
    f.add_scatter(
        "mapc_targets",
        {
            "x": x,
            "y": y,
            "c": [filtered_hac, filtered_frac_aromatic, filtered_num_rings, filtered_clogp, filtered_frac_csp3],
            "labels": labels,
        },
        point_scale=4,
        interactive=True,
        categorical=[False, False, False, False, False],
        colormap=['viridis', 'viridis', 'viridis', 'viridis', 'viridis'],
        series_title=['HAC', 'Fraction Aromatic Atoms', 'Number of Rings', 'clogP', 'Fraction Csp3'],
        has_legend=True,
    )

    # Add tree
    f.add_tree("mapc_targets_tree", {"from": s, "to": t}, point_helper="mapc_targets", color="#222222")
    
    # Plot
    f.plot('mapc_targets', template='smiles')


def main():
    
    data = pd.read_csv('data/230M/cluster_representatives.csv')

    # Extract PCA coordinates (fingerprints)
    pca_columns = ['PCA_1', 'PCA_2', 'PCA_3'] 
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

    # Set up the k-NN model for finding nearest neighbors in 5D space
    knn = 20
    print('Starting KNN')
    knn_search = NearestNeighbors(n_neighbors=knn, metric='manhattan', n_jobs=-1)  # Using 'manhattan' as the metric
    knn_search.fit(fingerprints)


    print('Creating layout')
    # Step 5: Create the Tmap layout configuration
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 20
    cfg.mmm_repeats = 2
    cfg.sl_extra_scaling_steps = 5
    cfg.k = 20
    cfg.sl_scaling_type = tm.RelativeToAvgLength

    print('Generating TMAP layout')
    # Step 6: Generate the Tmap layout using the edge list
    x_, y_, s, t, gp = tm.layout_from_edge_list(len(fingerprints), edge_list, cfg)
    
    plot_faerun(x_, y_, s, t, data)


if __name__ == "__main__":
    main()
