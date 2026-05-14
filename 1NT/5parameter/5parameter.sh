#!/bin/bash -l
#SBATCH -J 5param_all
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/5parameter/logs/5parameter_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/5parameter/logs/5parameter_%j.err

set -euo pipefail

PYBIN="${PYBIN:-$HOME/.conda/envs/work1/bin/python}"
RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/1NT/5parameter"

mkdir -p "$BASE/logs" "$BASE/data"

echo "[5parameter] full pipeline"
echo "[5parameter] partition: GPU-A40"
echo "[5parameter] cpus-per-task: 1"
echo "[5parameter] python: $PYBIN"
echo "[5parameter] rscript: $RBIN"

"$RBIN" "$BASE/prepare.R"
"$PYBIN" "$BASE/restart/cal_restart_mossn.py"
"$RBIN" "$BASE/restart/auc_restart.R"
"$RBIN" "$BASE/restart/plot_restart.R"
"$PYBIN" "$BASE/seed_threshold/cal_seed_threshold_mossn.py"
"$RBIN" "$BASE/seed_threshold/auc_seed_threshold.R"
"$RBIN" "$BASE/seed_threshold/plot_seed_threshold.R"
"$PYBIN" "$BASE/gamma/cal_gamma_mossn.py"
"$RBIN" "$BASE/gamma/auc_gamma.R"
"$RBIN" "$BASE/gamma/plot_gamma.R"
