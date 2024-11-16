#!/bin/bash

# Define paths
DATA_DIR="/mnt/samsung_2tb/mixed_data"
LOG_DIR="./logs"
SCRIPT_PATH="scripts/cluster_pipeline.py"

# Create the log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Log the total start time
echo "Total start time: $(date)" | tee -a "$LOG_DIR/total_time.log"
TOTAL_START=$(date +%s)

# Loop through files from 0 to 9 and run them sequentially
for i in {1..9}; do
    DATA_FILE="$DATA_DIR/output_file_$i.cxsmiles"
    OUTPUT_DIR="$DATA_DIR/output_file_$i"
    LOG_FILE="$LOG_DIR/output_file_$i.log"

    # Log the start time of the current run
    echo "Starting run for output_file_$i at $(date)" | tee -a "$LOG_FILE"
    START_TIME=$(date +%s)

    # Execute the command without background execution (&)
    nohup python "$SCRIPT_PATH" --data-file "$DATA_FILE" --output-dir "$OUTPUT_DIR" --chunksize 10000000 > "$LOG_FILE" 2>&1

    # Log the end time of the current run
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "Finished run for output_file_$i at $(date)" | tee -a "$LOG_FILE"
    echo "Duration for output_file_$i: $DURATION seconds" | tee -a "$LOG_FILE"

    # Log to the total time log
    echo "Run for output_file_$i completed in $DURATION seconds" | tee -a "$LOG_DIR/total_time.log"
done

# Log the total end time
TOTAL_END=$(date +%s)
TOTAL_DURATION=$((TOTAL_END - TOTAL_START))
echo "Total end time: $(date)" | tee -a "$LOG_DIR/total_time.log"
echo "Total duration for all runs: $TOTAL_DURATION seconds" | tee -a "$LOG_DIR/total_time.log"

# Notify the user
echo "All runs completed successfully."
