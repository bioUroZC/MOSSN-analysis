#!/bin/bash
#SBATCH -J PyC             # Job name
#SBATCH -N 1                        # Number of nodes
#SBATCH --ntasks-per-node=1        # One task per node
#SBATCH --mem=40G                  # Total memory
#SBATCH --exclude=node04
#SBATCH -o PyC.out            # Standard output file
#SBATCH -e PyC.err            # Standard error file

# Activate conda environment
source ~/.bashrc
conda activate work1

# Print environment info
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Using Python: $(which python)"

# Run the Python analysis script
python 2WRW.py
python 3SSN.py
python 4SWEET.py
