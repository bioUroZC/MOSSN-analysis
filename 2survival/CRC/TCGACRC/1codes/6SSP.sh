#!/bin/bash
#SBATCH -J 6S
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 6S.out
#SBATCH -e 6S.err

Rscript 6SSP.R
