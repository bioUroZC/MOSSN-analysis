#!/bin/bash -l
#SBATCH -J AUCcal
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o /proj/c.zihao/work1/1NT/2CPTAC/logs/4AUCcal_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/2CPTAC/logs/4AUCcal_%j.err

set -euo pipefail

DIR="/proj/c.zihao/work1/1NT/2CPTAC"
SCRIPT="$DIR/4AUCcal.R"

mkdir -p "$DIR/logs"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

[ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

echo "[start] $(date)"
echo "[info] node=$HOSTNAME cpus=$SLURM_CPUS_PER_TASK mem=16G"
echo "[info] running $SCRIPT"

srun --cpu-bind=cores Rscript "$SCRIPT"

echo "[done] $(date)"
