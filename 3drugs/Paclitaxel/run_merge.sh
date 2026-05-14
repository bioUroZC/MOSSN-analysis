#!/bin/bash -l
#SBATCH -J Paclitaxel_merge
#SBATCH -p GPU-A40
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/3drugs/Paclitaxel/logs/merge_%j.out
#SBATCH -e /proj/c.zihao/work1/3drugs/Paclitaxel/logs/merge_%j.err

set -euo pipefail

mkdir -p /proj/c.zihao/work1/3drugs/Paclitaxel/logs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[env] NODE=$HOSTNAME PART=${SLURM_JOB_PARTITION:-NA} CPUS=$SLURM_CPUS_PER_TASK"

srun --cpu-bind=cores Rscript /proj/c.zihao/work1/3drugs/Paclitaxel/2merge.R
