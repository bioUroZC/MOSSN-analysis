#!/bin/bash -l
#SBATCH -J GBM_lasso
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/2survival/logs/OS_GBM_%j.out
#SBATCH -e /proj/c.zihao/work1/2survival/logs/OS_GBM_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/2survival/GBM/2sur"
WORKDIR="$BASE"

echo "[GBM] prepare: $WORKDIR/1prepare.R"
"$RBIN" "$WORKDIR/1prepare.R"

echo "[GBM] lasso: $WORKDIR/2lassoSet.R"
"$RBIN" "$WORKDIR/2lassoSet.R"

echo "[GBM] lasso_rep: $WORKDIR/2lassoRep.R"
"$RBIN" "$WORKDIR/2lassoRep.R"

echo "[GBM] done"
