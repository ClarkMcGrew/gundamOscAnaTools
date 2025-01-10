
## Prepare Inputs for GUNDAM

This repository contains a script that prepares files in a GUNDAM-friendly format by copying the `cafTree` and adding new branches, including the converted dials into a `TClonesArray` of `TGraphs`.

### New Branches Added
The following branches are added to the `event_tree`:

- `POTScaledWeight`
- `POTweight`
- `Nonswap`
- `Nueswap`
- `Tauswap`
- `isNC`
- `isRHC`
- `Converted dials`

## How to Build and Run the Script

### Compiling the Script

To compile the `prepareGundamMCTree.cpp` script, use the following command:

```bash
g++ -o prepareGundamMCTree prepareGundamMCTree.cpp -I$(root-config --incdir) $(root-config --libs) -std=c++17
```

This command compiles the C++ code and links it against ROOT libraries. Make sure ROOT is properly installed and configured on your system.

### Running the Script Locally

Once the script is compiled, use the `run_jobs.sh` script to execute it. You may need to adjust the `useRHC` variable to `0` (for FHC mode) or `1` (for RHC mode) as needed.

Make the script executable and run it:

```bash
chmod +x run_jobs.sh
useRHC=1 # Adjust useRHC=1/0 as needed
splitCount=3
./run_jobs.sh $useRHC $splitCount
```

### Submitting Jobs to SLURM 
#### 1. Seawulf Cluster

To submit jobs to the Seawulf cluster using SLURM, use the following commands depending on the mode:

- For **FHC mode** (useRHC=0):

```bash
sbatch submit_slurm_FHC.sh
```

- For **RHC mode** (useRHC=1):

```bash
sbatch submit_slurm_RHC.sh
```

#### 2. NNhome
Submit the Job with a Dynamically Set Array:
```bash
useRHC=1 # Adjust useRHC=1/0 (RHC/FHC) as needed
splitCount=3 # Adjust job numbers as needed
sbatch --array=0-$((splitCount - 1)) nnhome_slurm.sh $useRHC $splitCount
```
