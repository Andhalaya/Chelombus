#!/usr/bin/env python3

import os
import sys
import xxhash

# Configuration
SEED = 42
NUM_OUTPUT_FILES = 10
OUTPUT_DIR = '/mnt/samsung_2tb/mixed_data'
OUTPUT_PREFIX = 'output_file_'
HASH_FUNC = xxhash.xxh64
BUFFER_SIZE = 16 * 1024 * 1024  # 16MB buffer

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def process_line(line, output_file_handles):
    molecule_smiles = line.strip()
    hash_value = HASH_FUNC(molecule_smiles.encode('utf-8'), seed=SEED).intdigest()
    file_index = hash_value % NUM_OUTPUT_FILES
    output_file_handles[file_index].write(f'{hash_value}\t{molecule_smiles}\n')

def process_input_file(input_file):
    output_file_paths = [os.path.join(OUTPUT_DIR, f'{OUTPUT_PREFIX}{i}.cxsmiles') for i in range(NUM_OUTPUT_FILES)]
    output_file_handles = [open(path, 'a', buffering=BUFFER_SIZE) for path in output_file_paths]

    with open(input_file, 'r', buffering=BUFFER_SIZE) as infile:
        for line in infile:
            process_line(line, output_file_handles)

    # Close output file handles
    for f in output_file_handles:
        f.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: sort_file.py input_file")
        sys.exit(1)
    input_file = sys.argv[1]
    process_input_file(input_file)
