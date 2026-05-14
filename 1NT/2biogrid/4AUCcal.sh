#!/bin/bash -l
#SBATCH -J biogrid_auc
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/2biogrid/logs/4AUCcal_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/2biogrid/logs/4AUCcal_%j.err

set -euo pipefail

DIR="/proj/c.zihao/work1/1NT/2biogrid"
SCRIPT="$DIR/4AUCcal.R"

mkdir -p "$DIR/logs"

cd "$DIR"

echo "[env] NODE=$HOSTNAME PART=$SLURM_JOB_PARTITION CPUS=$SLURM_CPUS_PER_TASK"
echo "[run] $(date '+%F %T') $SCRIPT"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpu-bind=cores Rscript "$SCRIPT"

echo "[done] $(date '+%F %T')"
