#!/bin/bash
#SBATCH -J 1download
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH -o 1download.out
#SBATCH -e 1download.err

Rscript 1download.R
