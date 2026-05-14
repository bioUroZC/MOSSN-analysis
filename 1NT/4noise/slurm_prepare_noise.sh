#!/bin/bash -l
#SBATCH -J 4noise_prep
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-3%4
#SBATCH -o /proj/c.zihao/work1/1NT/4noise/logs/prepare_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/4noise/logs/prepare_%A_%a.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/4noise"
LEVELS=(05 10 15 20)
LEVEL="${LEVELS[$SLURM_ARRAY_TASK_ID]}"

mkdir -p "$BASE/logs" "$BASE/$LEVEL/data"

echo "[4noise] prepare level=${LEVEL}"
echo "[4noise] rscript: $RBIN"

"$RBIN" "$BASE/$LEVEL/prepare.R"
