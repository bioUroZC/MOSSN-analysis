#!/bin/bash
#SBATCH -J download
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o RlC.out
#SBATCH -e RlC.err

Rscript 1download.R
