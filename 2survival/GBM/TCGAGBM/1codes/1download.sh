#!/bin/bash
#SBATCH -J 1download
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 1download.out
#SBATCH -e 1download.err

Rscript 1download.R

