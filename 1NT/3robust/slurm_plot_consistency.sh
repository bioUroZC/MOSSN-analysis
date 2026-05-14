#!/bin/bash -l
#SBATCH -J 3robust_plot
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -o /proj/c.zihao/work1/1NT/3robust/logs/plot_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/3robust/logs/plot_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/3robust"

mkdir -p "$BASE/logs" "$BASE/plots"

echo "[3robust] plot summary"
echo "[3robust] rscript: $RBIN"

"$RBIN" "$BASE/plot_consistency.R"
