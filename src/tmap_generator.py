import os 
import pandas as pd
import numpy as np

from typing import Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from faerun import Faerun
from config import OUTPUT_FILE_PATH
from fingerprint_calculator import FingerprintCalculator
from layout_computer import LayoutComputer, Plotter

class TmapConstructor: # Class to generate physicochemical properties from smiles 
    """
    properties: """
    def __init__(self, properties:tuple):
        pass

    def get_mol_properties(self, smiles):
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


class TmapGenerator:
    def __init__(
            self, 
            fingerprint_type: str = 'mhfp', 
            permutations: int=512, 
            method: str='lsh', 
            output_name: str = 'tmap', 
            labels: list =None, 
            vector_type: str = 'fingerprint', 
            cluster_id: Optional[str] = None):
        """
        type_vector: Select whether you want the fingerprints or the PCA coordinates to be used in the first TMAP
        cluster_id (Optional) Find the cluster_id which TMAP we will do. It has to be in format PCA1_PCA2_PCA3. 
        e.g. 0_12_10. Right now it finds the csv file with this label. In the future it will retrieve it from the database
        """

        self.fingerprint_type = fingerprint_type
        self.permutations = permutations
        self.method = method
        self.output_name = output_name
        self.labels = labels
        self.vector_type= vector_type 
        self._cluster_id = cluster_id
        # TODO: Think how to pass the correct dataframe
        self._dataframe = None
        
        # Initialize helper classes
        self.fingerprint_calculator = FingerprintCalculator(self._dataframe['smiles'], self.fingerprint_type, permutations=self.permutations)
        self.layout_computer = LayoutComputer(self.method, config=None)
        self.plotter = Plotter(self.output_name, self.labels)



        self.representatives_dataframe_file_path = os.path.join(OUTPUT_FILE_PATH, 'cluster_representatives.csv')
        self._representatives_dataframe = pd.read_csv(self.representatives_dataframe_file_path)
        

    def _get_fingerprint_vectors(self, encoder:str):
        """
        The vectors used for the TMAP layout will be the fingerprints calcualted from the SMILES
        """
        fp_calculator = FingerprintCalculator()  


    def _get_pca_vectors(self):
        pass


    def generate_representatives_tmap(self):
        """
        The vectors used for the TMAP layout will be the PCA values calculated previously
        To be used in the first TMAP (i.e. TMAPs of clusters using representatives as the nodes)
        """
        if self.type_vector == 'fingerprint':
            self.get_fingerprint_vectors

        elif self.type_vector == 'coordinates':
            self.get_pca_vectors

    def plot_farun(self):
        f = Faerun(view="front", 
                    coords=False, 
                    title= "Representatives Cluster", 
                    clear_color="#FFFFFF")
    def generate_cluster_tmap(self):
        """
        Generate the TMAP 'on the fly' based on the cluster_id given. It will look for the csv file for now but the idea is to retrieve from database
        """ 