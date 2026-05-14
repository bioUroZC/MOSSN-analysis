#!/bin/bash -l
#SBATCH -J merge_auc
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/2string/logs/merge_auc_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/2string/logs/merge_auc_%j.err

set -euo pipefail

DIR="/proj/c.zihao/work1/1NT/2string"
mkdir -p "$DIR/logs"
cd "$DIR"

echo "[env] NODE=${SLURMD_NODENAME:-unknown} PART=${SLURM_JOB_PARTITION:-unknown} CPUS=${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[run] $(date '+%F %T') 3merge.R"
srun --cpu-bind=cores Rscript "$DIR/3merge.R"
echo "[done] $(date '+%F %T') 3merge.R"

echo "[run] $(date '+%F %T') 4AUCcal.R"
srun --cpu-bind=cores Rscript "$DIR/4AUCcal.R"
echo "[done] $(date '+%F %T') 4AUCcal.R"
