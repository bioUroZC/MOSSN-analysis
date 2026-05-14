#!/bin/bash
#SBATCH -J 5l
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o 5l.out
#SBATCH -e 5l.err

Rscript 5lion.R
