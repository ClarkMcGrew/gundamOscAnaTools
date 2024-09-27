#!/bin/bash

# Load any required modules (if necessary)
# On SeaWulf
# module load root
# On NNhome machine
source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh

# Set useRHC to 0 (for FHC) or 1 (for RHC)
useRHC=${1:-0}  # Takes the first argument as useRHC, defaults to 0 if not provided

# Define the total number of splits (jobs) for manual loop
splitCount=${2:-1}  # Default splitCount is 1 if not provided

# Loop through the splits manually
for jobID in $(seq 0 $((splitCount - 1))); do
    echo "Running job with split ID $jobID"
    # Run the job in the background
    ./prepareGundamMCTree $useRHC $splitCount $jobID &
done

# Wait for all background jobs to finish
wait
