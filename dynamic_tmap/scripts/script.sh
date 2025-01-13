#!/usr/bin/env bash

# Number of parallel jobs, adjust based on your CPU cores
NPROC=4

# Run each line of cluster.txt in parallel
cat cluster.txt | parallel -j "$NPROC" 'time python tmap_pipeline.py --l {}'
