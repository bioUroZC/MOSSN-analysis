#!/bin/bash -l
#SBATCH -J pipe
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --array=0-4
#SBATCH -o logs/MOSSN_%A_%a.out
#SBATCH -e logs/MOSSN_%A_%a.err

set -euo pipefail

PYBIN="$HOME/.conda/envs/work1/bin/python"
DIR="/proj/c.zihao/work1/1NT/2intact"

export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "[env] NODE=$HOSTNAME PART=$SLURM_JOB_PARTITION CPUS=$SLURM_CPUS_PER_TASK"
"$PYBIN" -V

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

SCRIPTS=(
    "$DIR/cal_MOSSN_noPrior.py"
    "$DIR/cal_MOSSN_noSeed.py"
    "$DIR/cal_MOSSN_noRWR.py"
    "$DIR/cal_MOSSN_uniform.py"
    "$DIR/cal_MOSSN_noCorr.py"
)

SCRIPT=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
[ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

echo "[task ${SLURM_ARRAY_TASK_ID}] Running: $SCRIPT"

case "${SCRIPT##*.}" in
    py) srun --cpu-bind=cores "$PYBIN" "$SCRIPT" ;;
    R)  srun --cpu-bind=cores Rscript "$SCRIPT" ;;
    *)  echo "[ERR] Unknown extension: $SCRIPT"; exit 2 ;;
esac
