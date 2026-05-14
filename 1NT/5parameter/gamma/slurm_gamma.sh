#!/bin/bash -l
#SBATCH -J 5param_gam
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/5parameter/logs/gamma_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/5parameter/logs/gamma_%j.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/5parameter"

mkdir -p "$BASE/logs"

echo "[5parameter] gamma sensitivity"
echo "[5parameter] python: $PYBIN"
echo "[5parameter] rscript: $RBIN"

"$PYBIN" "$BASE/gamma/cal_gamma_mossn.py"
"$RBIN" "$BASE/gamma/auc_gamma.R"
"$RBIN" "$BASE/gamma/plot_gamma.R"
