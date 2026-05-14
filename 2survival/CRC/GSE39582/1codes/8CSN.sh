#!/bin/bash
#SBATCH -J 8C
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 8C.out
#SBATCH -e 8C.err

Rscript 8CSN.R
