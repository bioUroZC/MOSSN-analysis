#!/bin/bash -l
#SBATCH -J 4noise_mc
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --mem=24G
#SBATCH --array=0-3%4
#SBATCH -o /proj/c.zihao/work1/1NT/4noise/logs/merge_consistency_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/4noise/logs/merge_consistency_%A_%a.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/4noise"
LEVELS=(05 10 15 20)

LEVEL="${LEVELS[$SLURM_ARRAY_TASK_ID]}"
WORKDIR="$BASE/$LEVEL"

echo "[4noise] level=${LEVEL}"
echo "[4noise] merge: $WORKDIR/merge.R"
"$RBIN" "$WORKDIR/merge.R"

echo "[4noise] consistency: $WORKDIR/consistency.R"
"$RBIN" "$WORKDIR/consistency.R"

echo "[4noise] done: level=${LEVEL}"
