#!/bin/bash
#SBATCH -J RlC
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --exclude=node04
#SBATCH -o RlC.out
#SBATCH -e RlC.err

Rscript 5lion.R
Rscript 6SSP.R
Rscript 7iENA.R
Rscript 8CSN.R
