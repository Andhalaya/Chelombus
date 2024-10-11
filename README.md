# Chelombus

## Project Overview
Objective: Process a molecular dataset to compute fingerprints, reduce dimensionality using PCA, and output 3D coordinates along with SMILES strings and any available features.

Requirements:

    Adaptable to datasets with varying features or just SMILES.
    Output should include:
        3D coordinates for each molecule.
        Original SMILES strings.
        Available features from the dataset



### Structure

molecule-pca-visualization/
├── data/
│   └── sample_dataset.csv         # Place your datasets here
├── src/
│   ├── __init__.py
│   ├── data_handler.py            # Contains DataHandler class
│   ├── fingerprint_calculator.py  # Contains FingerprintCalculator class
│   ├── dimensionality_reducer.py  # Contains DimensionalityReducer class
│   ├── output_generator.py        # Contains OutputGenerator class
│   └── utils/
│       ├── __init__.py
│       └── helper_functions.py    # Any additional helper functions
├── notebooks/
│   └── exploratory_analysis.ipynb # Jupyter notebooks for testing and analysis
├── tests/
│   ├── __init__.py
│   ├── test_data_handler.py       # Unit tests for DataHandler
│   ├── test_fingerprint_calculator.py
│   ├── test_dimensionality_reducer.py
│   └── test_output_generator.py
├── scripts/
│   └── run_pipeline.py            # Script to execute the entire pipeline
├── requirements.txt               # Python dependencies
├── README.md                      # Project description and instructions
├── .gitignore                     # Files and folders to ignore in Git
└── LICENSE                        # License information
