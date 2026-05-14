#!/bin/bash -l
#SBATCH -J 4noise_plot
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -o /proj/c.zihao/work1/1NT/4noise/logs/plot_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/4noise/logs/plot_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/4noise"

mkdir -p "$BASE/logs" "$BASE/plots"

echo "[4noise] plot summary"
echo "[4noise] rscript: $RBIN"

"$RBIN" "$BASE/plot_consistency.R"
