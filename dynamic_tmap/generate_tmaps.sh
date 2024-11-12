#!/bin/bash

# Iterate over the first value (0 to 100)
for i in {0..100}; do
  # Iterate over the second value (0 to 100)
  for j in {0..100}; do
    # Iterate over the third value (0 to 50)
    for k in {0..50}; do
      label="${i}_${j}_${k}"
      output_dir="shared_volume/${i}/${j}/"

      # Create the output directory if it doesn't exist
      mkdir -p "$output_dir"

      # Run the pipeline script
      python scripts/tmap_pipeline.py --label "$label" --output_dir "$output_dir"
    done
  done
done
