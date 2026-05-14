#!/bin/bash
#SBATCH -J SSN            # Job name
#SBATCH -N 1                        # Number of nodes
#SBATCH --ntasks-per-node=1        # One task per node
#SBATCH --mem=40G                  # Total memory
#SBATCH -o SSN.out            # Standard output file
#SBATCH -e SSN.err            # Standard error file

# Activate conda environment
source ~/.bashrc
conda activate work1

# Print environment info
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Using Python: $(which python)"

# Run the Python analysis script
python 3SSN.py
