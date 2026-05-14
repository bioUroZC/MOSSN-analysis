#!/bin/bash -l
#SBATCH -J 3robust_mc
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --mem=24G
#SBATCH --array=0-3%4
#SBATCH -o /proj/c.zihao/work1/1NT/3robust/logs/merge_consistency_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/3robust/logs/merge_consistency_%A_%a.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/3robust"
LEVELS=(10 30 50 70)

LEVEL="${LEVELS[$SLURM_ARRAY_TASK_ID]}"
WORKDIR="$BASE/$LEVEL"

echo "[3robust] level=${LEVEL}"
echo "[3robust] merge: $WORKDIR/merge.R"
"$RBIN" "$WORKDIR/merge.R"

echo "[3robust] consistency: $WORKDIR/consistency.R"
"$RBIN" "$WORKDIR/consistency.R"

echo "[3robust] done: level=${LEVEL}"
