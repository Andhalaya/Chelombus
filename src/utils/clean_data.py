import re
import pandas as pd

def clean_smiles_file(input_file, output_file, is_csv):
    """
    Cleans the SMILES data in the input file by removing extraneous information such as '|&1:18,20|'.
    The cleaned data is saved to a new output file.
    
    :param input_file: Path to the input file (CSV or TXT/CXSMILES)
    :param output_file: Path to the output file where cleaned data will be saved
    :param is_csv: Whether the input file is a CSV (True) or a TXT/CXSMILES file (False)
    """
    def clean_smiles(smiles):
        """
        Removes any extraneous information from the SMILES string.
        Specifically, it removes patterns like '|&1:18,20|'.
        
        :param smiles: The raw SMILES string
        :return: The cleaned SMILES string
        """
        return re.sub(r'\|\&\d+:\d+(,\d+)*\|', '', smiles).strip()

    if is_csv:
        # Read CSV file
        data = pd.read_csv(input_file)

        # Clean the SMILES column
        smiles_column = data.filter(like='smiles').columns[0]  # Identify the SMILES column
        data[smiles_column] = data[smiles_column].apply(clean_smiles)  # Clean the SMILES column

        # Save the cleaned data to a new CSV file
        data.to_csv(output_file, index=False)
    else:
        # Read TXT or CXSMILES file line by line
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # Split line by tab
                parts = line.split('\t')
                if len(parts) >= 18:  # Ensure there are enough columns (SMILES + 17 features)
                    smiles = parts[0]
                    # Clean the SMILES string
                    cleaned_smiles = clean_smiles(smiles)
                    # Reconstruct the line with the cleaned SMILES
                    outfile.write(cleaned_smiles + '\n')

    print(f"Cleaning complete. Cleaned data saved to {output_file}.")


if __name__ == "__main__":
    # Example usage
    input_file = "data/650M.cxsmiles"  # Path to the original input file
    output_file =  "data/cleaned_650M.cxsmiles" # Path to the new file
    is_csv = False  # Set to True if input is a CSV, False for TXT/CXSMILES

    clean_smiles_file(input_file, output_file, is_csv)
