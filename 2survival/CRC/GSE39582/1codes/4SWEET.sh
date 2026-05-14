#!/bin/bash
#SBATCH -J 4S            # Job name
#SBATCH -N 1                        # Number of nodes
#SBATCH --ntasks-per-node=1        # One task per node
#SBATCH --mem=40G                  # Total memory
#SBATCH -o 4S.out            # Standard output file
#SBATCH -e 4S.err            # Standard error file

# Activate conda environment
source ~/.bashrc
conda activate work1

# Print environment info
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Using Python: $(which python)"

# Run the Python analysis script
python 4SWEET.py
