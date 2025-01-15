#!/usr/bin/env bash

# Set the number of parallel jobs based on Machine 1's cores
NPROC=6 # Adjust this based on the machine's optimal use

# Determine total lines and halfway point
TOTAL=$(wc -l < remaining_numbers.txt)
HALF=$(( (TOTAL * 62) / 100 ))
# Process the first half of the file using GNU Parallel
head -n "$HALF" remaining_numbers.txt | parallel -j "$NPROC" 'time python tmap_pipeline.py --l {}'
