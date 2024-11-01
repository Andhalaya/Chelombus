import os 
import pandas as pd
import numpy as np

from typing import Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from faerun import Faerun
from config import OUTPUT_FILE_PATH, TMAP_NODE_SIZE, TMAP_K, TMAP_POINT_SCALE
from src.fingerprint_calculator import FingerprintCalculator
from src.layout_computer import LayoutComputer, Plotter
import tmap as tm
import logging
 # Class to generate physicochemical properties from smiles 
class TmapConstructor:
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def _calculate_threshold(self, data):
        """ 
        Typically we could have very different values in the molecular properties (e.g. number of rings on a very diverse dataframe)
        which leads to lose of information due to outliers having extreme color values and the rest falling into the same range. 
        This function calculates a threshold using IQR method to cut the outliers value based on percentiles.
        """
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1 
        threshold = q3 + 1.5*iqr
        return threshold

    def _mol_properties_from_smiles(self, smiles: str) -> tuple:
        """ Get molecular properties from a single SMILES string"""
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
   
    def mol_properties_from_df(self): # -> list[list[float]]:
        self.dataframe[['hac', 'frac_aromatic', 'num_rings', 'clogp', 'frac_csp3']]  = self.dataframe['smiles'].apply(
            self._mol_properties_from_smiles
        ).apply(pd.Series)

        # Drop rows with any None or NaN values in the property columns
        df_clean = self.dataframe.dropna(subset=['hac', 'frac_aromatic', 'num_rings', 'clogp', 'frac_csp3'])
    
        # Calculate thresholds for each property using the clean DataFrame
        hac_threshold = self._calculate_threshold(df_clean['hac'])
        frac_threshold= self._calculate_threshold(df_clean['frac_aromatic'])
        rings_threshold=self._calculate_threshold(df_clean['num_rings'])
        clogp_threshold=self._calculate_threshold(df_clean['clogp'])
        csp3_threshold =self._calculate_threshold(df_clean['frac_csp3'])
    
        # Filter the DataFrame based on the thresholds
        filtered_df = df_clean[
            (df_clean['hac'] <= hac_threshold) &
            (df_clean['frac_aromatic'] <= frac_threshold) &
            (df_clean['num_rings'] <= rings_threshold) &
            (df_clean['clogp'] <= clogp_threshold) &
            (df_clean['frac_csp3'] <= csp3_threshold)
        ]
    
        # Extract filtered properties as lists
        filtered_hac = filtered_df['hac'].tolist()
        filtered_frac_aromatic = filtered_df['frac_aromatic'].tolist()
        filtered_num_rings = filtered_df['num_rings'].tolist()
        filtered_clogp = filtered_df['clogp'].tolist()
        filtered_frac_csp3 = filtered_df['frac_csp3'].tolist()
    
        # Return the list of lists
        return [filtered_hac, filtered_frac_aromatic, filtered_num_rings, filtered_clogp, filtered_frac_csp3]
    
 
class TmapGenerator:
    def __init__(
            self,
            dataframe: pd.DataFrame,  
            fingerprint_type: str = 'mhfp', 
            permutations: Optional[int]=512, 
            method: Optional[str]='lsh', 
            output_name: str = 'tmap', 
            categ_cols: list = None): 
        """
        param: fingerprint_type : Type of molecular fingerprint to be used on the TMAP. Options: {'mhfp', 'mqn', 'morgan', 'mapc'
        param: permutations: On MHFP number of permutations to be used in the MinHash

        TODO: Is KNN necessary? 
        param: method: Type of method to be used in TMAP. It can be either LSH (normal TMAP) or KNN? We could just create the TMAP always with LSH....
        
        param: output_name: name for the TMAP files. In case of the dynamic TMAP should inherit name from the cluster_id 
        param: categ_cols: List with the column names in the dataframe to be included as labels in the TMAP. These are typically your categorical columns
        e.g. 'Target_type', 'Organism' ...  
        """


        self.fingerprint_type = fingerprint_type
        self.permutations = permutations
        self.method = method
        self.output_name = output_name
        self.categ_cols = categ_cols 
        # TODO: Think how to pass the correct dataframe
        self.dataframe = dataframe 
        
        # Initialize helper classes
        self.fingerprint_calculator = FingerprintCalculator(self.dataframe['smiles'], self.fingerprint_type, permutations=self.permutations)
        self.layout_computer = LayoutComputer(self.method, config=None)
        self.plotter = Plotter(self.output_name, self.categ_cols)
        self.tmap_constructor = TmapConstructor(self.dataframe)

        #TODO: Is this necessary? I could just create a new instance of the class and treat as just any other TMAP passing 'cluster_id' as label'
        self.representatives_dataframe_file_path = os.path.join(OUTPUT_FILE_PATH, 'cluster_representatives.csv')
        self._representatives_dataframe = pd.read_csv(self.representatives_dataframe_file_path)
        

    def _get_fingerprint_vectors(self):
        """
        The vectors used for the TMAP layout will be the fingerprints calcualted from the SMILES
        """
        pass


    def _get_pca_vectors(self):
        """
        The vectors used for the TMAP layout will be the PCA Components. Typically to be used for the representatives cluster TMAP
        """
        #TODO: Get 3D Coordinates from self.representative_dataframe
        # Then use those as vectors for the TMAP. Not sure if worth it. Check with JL

    def generate_representatives_tmap(self, type_vector: str='fingerprint'):
        """
        Hierarchical TMAP of the clusters. Each point/node in the TMAP is a representative molecule of the cluster. 
        The vectors used for the TMAP layout will be either the 'fingerprint' of the representative molecule 
        or the PCA values (3D Coordinates). Default is 'fingerprint' 
        """
        if type_vector == 'fingerprint':
            self._get_fingerprint_vectors()
        elif type_vector == 'coordinates':
            self._get_pca_vectors()
        self.label = 'cluster_id'

    def construct_lsh_forest(self) -> None:
        fingerprints = self.fingerprint_calculator.calculate_fingerprints() 
        tm_fingerprints  = [tm.VectorUint(fp) for fp in fingerprints] #TMAP requires fingerprints to be passed as VectorUint

        lf = tm.LSHForest(512, 128, store=True)
        lf.batch_add(tm_fingerprints)
        lf.index()

        cfg = tm.LayoutConfiguration()
        cfg.node_size = TMAP_NODE_SIZE
        cfg.mmm_repeats = 2
        cfg.sl_extra_scaling_steps = 10
        cfg.k = TMAP_K
        cfg.sl_scaling_type = tm.RelativeToAvgLength
        
        logging.info("Layout")
        self.x, self.y, self.s, self.t, _ = tm.layout_from_lsh_forest(lf, cfg)

    def plot_faerun(self):
        self.construct_lsh_forest() 
        f = Faerun(view="front", 
                    coords=False, 
                    title= "Representatives Cluster", 
                    clear_color="#FFFFFF")

        # Create categories
        labels = []
        for i, row in self._dataframe.iterrows():
            label = '__'.join(str(row[col]) for col in self.labels)
            labels.append(label)
        
        properties = self.tmap_constructor.mol_properties_from_df() 

        # Categorical = [True] * categorical_columns + [False]*numerical_columns 
        numerical_col= [False]*5 # These are 5 by default. 5 molecular properties 
        categorical_col = [True]*len(self.categ_cols)
        bool_categorical = numerical_col + categorical_col

        colormap = ['tab10' if value else 'viridis' for value in bool_categorical]

        series_title = self.categ_cols + ['HAC', 'Fraction Aromatic Atoms', 'Number of Rings', 'clogP', 'Fraction Csp3']
        
        f.add_scatter(
            "mapc_targets", 
            {
                "x": self.x, 
                "y": self.y, 
                "c": properties, 
                "labels": labels, 
            }, 
        shader="smoothCircle",
        point_scale= TMAP_POINT_SCALE,
        max_point_size=20,
        interactive=True,
        legend_labels=[], # TODO: Get list of list of labels. This sould be something like [df[col] for col in self.categ_col]
        categorical= bool_categorical, 
        colormap= colormap, 
        series_title= series_title, 
        has_legend=True,
        )

        # Add tree
        f.add_tree("mapc_targets_tree", {"from": self.s, "to": self.t}, point_helper="mapc_targets", color="#222222")
        
        # Plot
        f.plot('mapc_targets', template='smiles')

    def generate_cluster_tmap(self):
        """
        Generate the TMAP 'on the fly' based on the cluster_id given. It will look for the csv file for now but the idea is to retrieve from database
        cluster_id (str int_int_int) Find the cluster_id which TMAP we will do. It has to be in format PCA1_PCA2_PCA3. 
        e.g. 0_12_10. Right now it finds the csv file with this label. In the future it will retrieve it from the database
        """
        pass