#!/bin/bash
#SBATCH -J 7i
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 7i.out
#SBATCH -e 7i.err

Rscript 7iENA.R
