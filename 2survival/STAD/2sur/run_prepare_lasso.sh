#!/bin/bash -l
#SBATCH -J STAD_lasso
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/2survival/logs/OS_STAD_%j.out
#SBATCH -e /proj/c.zihao/work1/2survival/logs/OS_STAD_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/2survival/STAD/2sur"
WORKDIR="$BASE"

echo "[STAD] prepare: $WORKDIR/1prepare.R"
"$RBIN" "$WORKDIR/1prepare.R"

echo "[STAD] lasso: $WORKDIR/2lassoSet.R"
"$RBIN" "$WORKDIR/2lassoSet.R"

echo "[STAD] lasso_rep: $WORKDIR/2lassoRep.R"
"$RBIN" "$WORKDIR/2lassoRep.R"

echo "[STAD] done"
