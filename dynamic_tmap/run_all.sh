#!/bin/bash

# Record the start time
start_time=$(date +%s)

# Path to the cluster_pipeline.py script
SCRIPT_DIR="$(dirname "$0")/scripts"
PYTHON_SCRIPT="$SCRIPT_DIR/cluster_pipeline.py"

# Path to the data directories
DATA_ROOT="/mnt/10tb_hdd/clickhouse_chelombus_data"

echo "DATA_ROOT: $DATA_ROOT"
echo "Directories found:"
ls "$DATA_ROOT"

# Check if the script has execute permissions
if [ ! -x "$0" ]; then
    echo "Script does not have execute permissions. Setting execute permissions."
    chmod +x "$0"
fi

# Verify that the Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Python script not found at $PYTHON_SCRIPT"
    exit 1
fi

# Loop over data directories
for dir in "$DATA_ROOT"/*/; do
    echo "Processing directory $dir"

    # Find the .cxsmiles file in the directory
    datafile=$(find "$dir" -maxdepth 1 -type f -name "*.cxsmiles")

    if [ -z "$datafile" ]; then
        echo "No .cxsmiles file found in $dir"
        continue
    fi

    # Check if datafile is readable
    if [ ! -r "$datafile" ]; then
        echo "Datafile $datafile is not readable. Setting read permissions."
        chmod +r "$datafile"
    fi

    # Check if the output directory is writable
    if [ ! -w "$dir" ]; then
        echo "Directory $dir is not writable. Setting write permissions."
        chmod +w "$dir"
    fi

    echo "Found datafile: $datafile"

    # Run the python script
    python "$PYTHON_SCRIPT" --data-file "$datafile" --output-dir "$dir"

done

# Record the end time
end_time=$(date +%s)

# Calculate the elapsed time in seconds
elapsed=$((end_time - start_time))

# Add 280.85 minutes of preprocessing time (converted to seconds)
additional_time_seconds=$(echo "280.85 * 60" | bc)
total_elapsed=$((elapsed + ${additional_time_seconds%.*}))

# Convert total elapsed time to hours, minutes, and seconds
hours=$((total_elapsed / 3600))
minutes=$(((total_elapsed % 3600) / 60))
seconds=$((total_elapsed % 60))

# Print the total time
echo "Total time including preprocessing: ${hours}h ${minutes}m ${seconds}s"
