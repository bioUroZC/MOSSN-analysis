#!/bin/bash -l
#SBATCH -J 5param_seed
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/5parameter/logs/seed_threshold_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/5parameter/logs/seed_threshold_%j.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/5parameter"

mkdir -p "$BASE/logs"

echo "[5parameter] seed-threshold sensitivity"
echo "[5parameter] python: $PYBIN"
echo "[5parameter] rscript: $RBIN"

"$PYBIN" "$BASE/seed_threshold/cal_seed_threshold_mossn.py"
"$RBIN" "$BASE/seed_threshold/auc_seed_threshold.R"
"$RBIN" "$BASE/seed_threshold/plot_seed_threshold.R"
