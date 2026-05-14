#!/bin/bash
#SBATCH -J 1Prepare
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 1Prepare.out
#SBATCH -e 1Prepare.err

Rscript 1Prepare.R
