#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --job-name=SplitMC
#SBATCH --output=slurm_files/slurm-%A_%a.out
#SBATCH --error=slurm_files/slurm-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --array=0-9    # 10 jobs (you can adjust this)
#SBATCH --time=02:00:00

module load root

# Set useRHC to 0 (for FHC) or 1 (for RHC)
useRHC=1  # Adjust accordingly

# Define the total number of splits (jobs) for parallel processing
splitCount=10  # You can change this value

# Call the program with useRHC, the total number of splits, and the current job ID
./prepareGundamMCTree $useRHC $splitCount $SLURM_ARRAY_TASK_ID
