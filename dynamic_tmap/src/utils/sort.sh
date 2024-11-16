#!/bin/bash

# Paths
INPUT_FILE='/mnt/10tb_hdd/uncleaned_data/Enamine_REAL_HAC_26_766M_CXSMILES.cxsmiles'
PROCESS_SCRIPT='./sort_file.py'

# Ensure the processing script is executable
chmod +x $PROCESS_SCRIPT

# Run the processing script
python3 $PROCESS_SCRIPT $INPUT_FILE
