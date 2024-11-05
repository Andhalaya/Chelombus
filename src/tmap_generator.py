import random 
import os 
from typing import List
import pandas as pd
import numpy as np

from typing import Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from faerun import Faerun
from config import OUTPUT_FILE_PATH, TMAP_NAME, TMAP_NODE_SIZE, TMAP_K, TMAP_POINT_SCALE
from src.fingerprint_calculator import FingerprintCalculator
from src.layout_computer import LayoutComputer, Plotter
import tmap as tm
import logging

 # Class to generate physicochemical properties from smiles 
class TmapConstructor:
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def _calculate_threshold(self, data) -> float:
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
        molecular_weight = rdMolDescriptors.ExactMolWt(mol)
        clogP = Descriptors.MolLogP(mol)
        fraction_Csp3 = Descriptors.FractionCSP3(mol)
        
        return (hac, fraction_aromatic_atoms, number_of_rings, clogP, fraction_Csp3)
   
    def mol_properties_from_df(self) -> List[List[float]]:
        self.dataframe[['hac', 'frac_aromatic', 'num_rings', 'clogp', 'frac_csp3']]  = self.dataframe['smiles'].apply(
            self._mol_properties_from_smiles
        ).apply(pd.Series)

        # Drop rows with any None or NaN values in the property columns
        # df_clean = self.dataframe.dropna(subset=['hac', 'frac_aromatic', 'num_rings', 'clogp', 'frac_csp3'])
    
        # Calculate thresholds for each property using the clean DataFrame
        hac_threshold = self._calculate_threshold(self.dataframe['hac'])
        frac_threshold= self._calculate_threshold(self.dataframe['frac_aromatic'])
        rings_threshold=self._calculate_threshold(self.dataframe['num_rings'])
        clogp_threshold=self._calculate_threshold(self.dataframe['clogp'])
        csp3_threshold =self._calculate_threshold(self.dataframe['frac_csp3'])
    
        # Filter the DataFrame based on the thresholds
        filtered_df = self.dataframe[
            (self.dataframe['hac'] <= hac_threshold) &
            (self.dataframe['frac_aromatic'] <= frac_threshold) &
            (self.dataframe['num_rings'] <= rings_threshold) &
            (self.dataframe['clogp'] <= clogp_threshold) &
            (self.dataframe['frac_csp3'] <= csp3_threshold)
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
            df_path: pd.DataFrame,  
            fingerprint_type: str = 'morgan', 
            permutations: Optional[int]= 512, 
            output_name: str = TMAP_NAME,
            fp_size: int = 1024, 
            categ_cols: Optional[List] = None
        ):
        """
        param: fingerprint_type : Type of molecular fingerprint to be used on the TMAP. Options: {'mhfp', 'mqn', 'morgan', 'mapc'
        param: permutations: On MHFP number of permutations to be used in the MinHash

        TODO: Is KNN necessary? 
        param: method: Type of method to be used in TMAP. It can be either LSH (normal TMAP) or KNN? We could just create the TMAP always with LSH....
        
        param: output_name: name for the TMAP files. In case of the dynamic TMAP should inherit name from the cluster_id 
        param: categ_cols: List with the column names in the dataframe to be included as labels in the TMAP. These are typically your categorical columns
        e.g. 'Target_type', 'Organism' ...  
        """
        self.df_path = df_path
        self.fingerprint_type = fingerprint_type
        self.permutations = permutations
        self.output_name = output_name
        self.fp_size = fp_size
        self.categ_cols = categ_cols 
        self.dataframe = pd.read_csv(df_path) 
        
        # Initialize helper classes
        self.fingerprint_calculator = FingerprintCalculator(self.dataframe['smiles'], self.fingerprint_type, permutations=self.permutations, fp_size=self.fp_size)
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

    def tmap_little(self): # Generates TMAP with minimum configuration, for testing purposes
        
        fingerprints = self.fingerprint_calculator.calculate_fingerprints()

        self.construct_lsh_forest(fingerprints)

        labels = []
        for i, row in self.dataframe.iterrows():
            if self.categ_cols != None:
                label = '__'.join(str(row[col]) for col in self.categ_cols)
                labels.append(row['smiles']+'__'+label)
            else:
                label.append(row['smiles'])

        # Plotting
        f = Faerun(
            view="front",
            coords=False,
            title="",
            clear_color="#FFFFFF",
        )

        f.add_scatter(
            TMAP_NAME+"_TMAP",
            {
                "x": self.x,
                "y": self.y,
                "c": np.arange(len(fingerprints)), # random values
                "labels":labels,
            },
            shader="smoothCircle",
            point_scale= TMAP_POINT_SCALE,
            max_point_size= 20,
            interactive=True,
            # legend_labels=[], # TODO: Get list of list of labels. This sould be something like [df[col] for col in self.categ_col]
            # categorical= bool_categorical, 
            # series_title= series_title, 
            has_legend=True,           
            colormap=['rainbow'],
            categorical=[False],
        )

        f.add_tree(TMAP_NAME+"_TMAP_tree", {"from": self.s, "to": self.t}, point_helper=TMAP_NAME+"_TMAP")
        f.plot(TMAP_NAME+"_TMAP", template='smiles')


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

    def construct_lsh_forest(self, fingerprints) -> None:
        tm_fingerprints  = [tm.VectorUint(fp) for fp in fingerprints] #TMAP requires fingerprints to be passed as VectorUint

        # LSH Indexing and coordinates generation
        logging.info("Plotting...")
        lf = tm.LSHForest(self.permutations)
        lf.batch_add(tm_fingerprints)
        lf.index()

        # Get the coordinates and Layout Configuration
        cfg = tm.LayoutConfiguration()
        cfg.node_size = TMAP_NODE_SIZE 
        cfg.mmm_repeats = 2
        cfg.sl_extra_scaling_steps = 10
        cfg.k = TMAP_K 
        cfg.sl_scaling_type = tm.RelativeToAvgLength
        self.x, self.y, self.s, self.t, _ = tm.layout_from_lsh_forest(lf, cfg)

    def plot_faerun(self, fingerprints):
        logging.info("Constructing LSH Forest...")
        self.construct_lsh_forest(fingerprints) 
        f = Faerun(view="front", 
                    coords=False, 
                    title= "", 
                    clear_color="#FFFFFF")
        
        def safe_create_categories(series):
            return Faerun.create_categories(series.fillna('Unknown').astype(str))
        
        
        # Create categories
        labels = []
        for i, row in self.dataframe.iterrows():
            if self.categ_cols != None:
                label = '__'.join(str(row[col]) for col in self.categ_cols)
                labels.append(row['smiles']+'__'+label)
            else:
                labels.append(row['smiles'])

        logging.info("Plotting...")
        properties = self.tmap_constructor.mol_properties_from_df() 

        # Categorical = [True] * categorical_columns + [False]*numerical_columns 
        numerical_col= [False]*5 # These are 5 by default. 5 molecular properties 
        categorical_col = [True]*len(self.categ_cols)
        bool_categorical = categorical_col + numerical_col  # List of booleans required to indicate if a label is categorical or numerical
        
        cluster_labels, cluster_data = safe_create_categories(self.dataframe['cluster_id'])

        colormap = ['tab10' if value else 'viridis' for value in bool_categorical]
        # properties.insert(0, cluster_data)
        series_title = self.categ_cols + ['HAC', 'Fraction Aromatic Atoms', 'Number of Rings', 'clogP', 'Fraction Csp3']
        c = random.sample(range(1,len(properties[2])*2), len(properties[2]))
        f.add_scatter(
            "mapc_targets", 
            {
                "x": self.x, 
                "y": self.y, 
                "c": c, 
                "labels": self.dataframe['smiles'], 
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
        
        f.add_tree("mhfp_tmap_node140_TMAP_tree", {"from": self.s, "to": self.t}, point_helper="mhfp_tmap_node140_TMAP")
        f.plot('mhfp_tmap_node140_TMAP', template='smiles')
        # Plot

    def generate_cluster_tmap(self):
        """
        Generate the TMAP 'on the fly' based on the cluster_id given. It will look for the csv file for now but the idea is to retrieve from database
        cluster_id (str int_int_int) Find the cluster_id which TMAP we will do. It has to be in format PCA1_PCA2_PCA3. 
        e.g. 0_12_10. Right now it finds the csv file with this label. In the future it will retrieve it from the database
        """
        pass