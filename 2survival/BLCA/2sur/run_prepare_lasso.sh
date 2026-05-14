#!/bin/bash -l
#SBATCH -J BLCA_lasso
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/2survival/logs/OS_BLCA_%j.out
#SBATCH -e /proj/c.zihao/work1/2survival/logs/OS_BLCA_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/2survival/BLCA/2sur"
WORKDIR="$BASE"

echo "[BLCA] prepare: $WORKDIR/1prepare.R"
"$RBIN" "$WORKDIR/1prepare.R"

echo "[BLCA] lasso: $WORKDIR/2lassoSet.R"
"$RBIN" "$WORKDIR/2lassoSet.R"

echo "[BLCA] lasso_rep: $WORKDIR/2lassoRep.R"
"$RBIN" "$WORKDIR/2lassoRep.R"

echo "[BLCA] done"
