#!/bin/bash -l
#SBATCH -J 5param_rst
#SBATCH -N 1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/5parameter/logs/restart_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/5parameter/logs/restart_%j.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/5parameter"

mkdir -p "$BASE/logs"

echo "[5parameter] restart sensitivity"
echo "[5parameter] python: $PYBIN"
echo "[5parameter] rscript: $RBIN"

"$PYBIN" "$BASE/restart/cal_restart_mossn.py"
"$RBIN" "$BASE/restart/auc_restart.R"
"$RBIN" "$BASE/restart/plot_restart.R"
