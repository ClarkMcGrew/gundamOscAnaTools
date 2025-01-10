#!/bin/bash
#SBATCH --job-name=SplitMC
#SBATCH --output=slurm_files/SplitMC-%A_%a.out
#SBATCH --error=slurm_files/SplitMC-%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Flynn.Y.Guo@stonybrook.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh

# Set useRHC to 0 (for FHC) or 1 (for RHC)
useRHC=${1:-0}  # Takes the first argument as useRHC, defaults to 0 if not provided

# Define the total number of splits (jobs) for parallel processing
splitCount=${2:-1}  # Default splitCount is 1 if not provided

# Call the program with useRHC, the total number of splits, and the current job ID
./prepareGundamMCTree $useRHC $splitCount $SLURM_ARRAY_TASK_ID
