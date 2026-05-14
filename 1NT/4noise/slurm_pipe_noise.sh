#!/bin/bash -l
#SBATCH -J 4noise_pipe
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --array=0-27%8
#SBATCH -o /proj/c.zihao/work1/1NT/4noise/logs/pipe_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/4noise/logs/pipe_%A_%a.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/4noise"
LEVELS=(05 10 15 20)

METHOD_SCRIPTS=(
  "cal_SWEET.py"
  "cal_SSN.py"
  "cal_LIONESS.R"
  "cal_MOSSN_noPrior.py"
  "cal_Patkar.py"
  "cal_PPIXpress.py"
  "cal_Proteinarium.py"
)

METHOD_COUNT=${#METHOD_SCRIPTS[@]}
level_idx=$((SLURM_ARRAY_TASK_ID / METHOD_COUNT))
script_idx=$((SLURM_ARRAY_TASK_ID % METHOD_COUNT))
LEVEL="${LEVELS[$level_idx]}"
DIR="$BASE/$LEVEL"
SCRIPT="$DIR/${METHOD_SCRIPTS[$script_idx]}"

export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

mkdir -p "$BASE/logs" "$DIR/logs"

echo "[env] level=$LEVEL node=$HOSTNAME part=$SLURM_JOB_PARTITION cpus=$SLURM_CPUS_PER_TASK"
echo "[env] python: $PYBIN"
echo "[env] rscript: $RBIN"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

[ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

echo "[task ${SLURM_ARRAY_TASK_ID}] Running: $SCRIPT"

case "${SCRIPT##*.}" in
  py) srun --cpu-bind=cores "$PYBIN" "$SCRIPT" ;;
  R)  srun --cpu-bind=cores "$RBIN" "$SCRIPT" ;;
  *)  echo "[ERR] Unknown extension: $SCRIPT"; exit 2 ;;
esac
