#!/bin/bash
#SBATCH -J run_all_features
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH -o run_all.out
#SBATCH -e run_all.err

Rscript 1wrwFeature.R > 1wrw.out 2> 1wrw.err
Rscript 2lionFeature.R > 2lion.out 2> 2lion.err
Rscript 3ssnFeature.R > 3ssn.out 2> 3ssn.err
Rscript 4ppixFeature.R > 4ppix.out 2> 4ppix.err
