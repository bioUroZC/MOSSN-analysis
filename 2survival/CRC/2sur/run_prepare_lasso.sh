#!/bin/bash -l
#SBATCH -J CRC_lasso
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/2survival/logs/OS_CRC_%j.out
#SBATCH -e /proj/c.zihao/work1/2survival/logs/OS_CRC_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/2survival/CRC/2sur"
WORKDIR="$BASE"

echo "[CRC] prepare: $WORKDIR/1prepare.R"
"$RBIN" "$WORKDIR/1prepare.R"

echo "[CRC] lasso: $WORKDIR/2lassoSet.R"
"$RBIN" "$WORKDIR/2lassoSet.R"

echo "[CRC] lasso_rep: $WORKDIR/2lassoRep.R"
"$RBIN" "$WORKDIR/2lassoRep.R"

echo "[CRC] done"
