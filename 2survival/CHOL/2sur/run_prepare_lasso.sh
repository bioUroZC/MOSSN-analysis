#!/bin/bash -l
#SBATCH -J CHOL_lasso
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/2survival/logs/OS_CHOL_%j.out
#SBATCH -e /proj/c.zihao/work1/2survival/logs/OS_CHOL_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/2survival/CHOL/2sur"
WORKDIR="$BASE"

echo "[CHOL] prepare: $WORKDIR/1prepare.R"
"$RBIN" "$WORKDIR/1prepare.R"

echo "[CHOL] lasso: $WORKDIR/2lassoSet.R"
"$RBIN" "$WORKDIR/2lassoSet.R"

echo "[CHOL] lasso_rep: $WORKDIR/2lassoRep.R"
"$RBIN" "$WORKDIR/2lassoRep.R"

echo "[CHOL] done"
