#!/bin/bash -l
#SBATCH -J 3robust_prep
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-3%4
#SBATCH -o /proj/c.zihao/work1/1NT/3robust/logs/prepare_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/3robust/logs/prepare_%A_%a.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/3robust"
LEVELS=(10 30 50 70)
LEVEL="${LEVELS[$SLURM_ARRAY_TASK_ID]}"

mkdir -p "$BASE/logs" "$BASE/$LEVEL/data"

echo "[3robust] prepare level=${LEVEL}"
echo "[3robust] rscript: $RBIN"

"$RBIN" "$BASE/$LEVEL/prepare.R"
