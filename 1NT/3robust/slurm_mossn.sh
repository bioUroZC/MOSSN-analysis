#!/bin/bash -l
#SBATCH -J 3robust_mossn_noprior
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --array=0-3%4
#SBATCH -o /proj/c.zihao/work1/1NT/3robust/logs/mossn_noprior_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/3robust/logs/mossn_noprior_%A_%a.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
BASE="/proj/c.zihao/work1/1NT/3robust"
LEVELS=(10 30 50 70)
SCRIPT_NAME="cal_MOSSN_noPrior.py"

LEVEL="${LEVELS[$SLURM_ARRAY_TASK_ID]}"
DIR="$BASE/$LEVEL"
SCRIPT="$DIR/$SCRIPT_NAME"

export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

mkdir -p "$BASE/logs" "$DIR/logs"

echo "[env] level=$LEVEL node=$HOSTNAME part=$SLURM_JOB_PARTITION cpus=$SLURM_CPUS_PER_TASK"
echo "[env] python: $PYBIN"
echo "[task ${SLURM_ARRAY_TASK_ID}] Running: $SCRIPT"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

[ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

srun --cpu-bind=cores "$PYBIN" "$SCRIPT"
