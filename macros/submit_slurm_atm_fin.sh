#!/bin/bash
#SBATCH -p extended-40core-shared
#SBATCH --job-name=SplitMC
#SBATCH --output=slurm_files/slurm-%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=0-35    # 36 jobs (you can adjust this)
#SBATCH --time=02:00:00

module load root

# Call the program with useRHC, the total number of splits, and the current job ID
./fileConvertor_ATM $SLURM_ARRAY_TASK_ID
