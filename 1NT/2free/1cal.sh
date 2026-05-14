#!/bin/bash -l
#SBATCH -J scalefree_string
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/2free/logs/scalefree_string_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/2free/logs/scalefree_string_%j.err

set -euo pipefail
mkdir -p /proj/c.zihao/work1/1NT/2free/logs
cd /proj/c.zihao/work1/1NT/2free

echo "[env] NODE=$HOSTNAME PART=$SLURM_JOB_PARTITION CPUS=$SLURM_CPUS_PER_TASK"
echo "[run] $(date '+%F %T') STRING"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun --cpu-bind=cores Rscript 1cal.R
echo "[done] $(date '+%F %T')"
