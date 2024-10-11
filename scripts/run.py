from src.data_handler import DataHandler
from src.fingerprint_calculator import FingerprintCalculator
from src.dimensionality_reducer import DimensionalityReducer
from src.output_generator import OutputGenerator


def main():
    # Initialize classes
    data_handler = DataHandler()
    fp_calculator = FingerprintCalculator()
    reducer = DimensionalityReducer()
    output_gen = OutputGenerator()

    # Step 1: Load and preprocess data
    data_handler.load_data('data/sample_dataset.csv')
    smiles_list = data_handler.extract_smiles()
    features = data_handler.extract_features()

    # Step 2: Calculate fingerprints
    fingerprints = fp_calculator.calculate_fingerprints(smiles_list)

    # Step 3: Dimensionality Reduction
    coordinates = reducer.fit_transform(fingerprints)

    # Step 4: Compile and export data
    output_gen.compile_output(coordinates, smiles_list, features)
    output_gen.export_data('data/processed_data.csv')

if __name__ == '__main__':
    main()
