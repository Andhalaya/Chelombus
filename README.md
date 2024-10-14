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

```bash
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
```


### PCA 

Initializing StandardScaler and IncrementalPCA:

StandardScaler() is initialized without parameters, as with_mean and with_std are True by default.
IncrementalPCA(n_components=n_components) is initialized with the desired number of components.
Step 1: Fitting the Scaler:

scaler.partial_fit(fingerprints_array) is called for each chunk to update the running mean and variance.
This step effectively computes the global mean and standard deviation across all chunks.
Step 2: Fitting Incremental PCA:

After the scaler has been fitted to the entire dataset, we standardize each chunk using scaler.transform(fingerprints_array).
We then reducer.partial_fit(fingerprints_std) to incrementally fit the PCA model on the standardized data.
Step 3: Transforming and Saving Data:

We standardize and transform each chunk using the fitted scaler and PCA model.
The reduced data is then saved to disk using pickle.
Important Considerations
Order of Operations: It's crucial to first fit the scaler to the entire dataset (or an approximation using chunks) before standardizing and fitting PCA. This ensures that the standardization parameters are consistent across all data.

Multiple Passes Over Data:

First Pass: Fit the StandardScaler to compute the mean and variance.
Second Pass: Fit the IncrementalPCA using the standardized data.
Third Pass: Transform the data using the fitted PCA model and save the results.
Unfortunately, this requires multiple passes over the data. If reading the data multiple times is not feasible, you might need to consider approximations or advanced techniques.

Memory Management: Even though we are processing data in chunks, ensure that each chunk is small enough to fit in memory, considering the overhead of additional data structures.

#### Alternative Approach: Combining Standardization and PCA in One Pass
If making multiple passes over the data is impractical, also we can consider the following approach:

Estimate Scaling Parameters with a Subset:
    Use a representative subset of your data to compute approximate mean and standard deviation.
    Fit the StandardScaler on this subset.
    Proceed with Single-Pass Processing:

For each chunk:
Standardize using the approximate scaler.
Partial fit the IncrementalPCA.
Optionally, transform and save the reduced data


## TO DO
- Add arguments parser to num_dimensions to reduce to 
- Add support for txt and other file typess
- Think of better output type than pd.DataFrame: options Parquet, HDF5, 
- Mkdir output in output_generator.py
- Add directory to save log output s