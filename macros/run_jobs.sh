#!/bin/bash

# Load any required modules (if necessary)
module load root

# Set useRHC to 0 (for FHC) or 1 (for RHC)
useRHC=1  # Adjust this value as needed

# Define the total number of splits (jobs) for manual loop
splitCount=10  # You can change this value

# Loop through the splits manually
for jobID in $(seq 0 $((splitCount - 1))); do
    echo "Running job with split ID $jobID"
    # Run the job in the background
    ./prepareGundamMCTree $useRHC $splitCount $jobID &
done

# Wait for all background jobs to finish
wait

