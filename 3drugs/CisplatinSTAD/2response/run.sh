#!/bin/bash -l
#SBATCH -J CisplatinSTAD_response
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o /proj/c.zihao/work1/3drugs/logs/CisplatinSTAD_response_%j.out
#SBATCH -e /proj/c.zihao/work1/3drugs/logs/CisplatinSTAD_response_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/3drugs/CisplatinSTAD/2response"

echo "[CisplatinSTAD] prepare: $BASE/1prepare.R"
"$RBIN" "$BASE/1prepare.R"

echo "[CisplatinSTAD] rf_dataset: $BASE/2rf_dataset.R"
"$RBIN" "$BASE/2rf_dataset.R"

echo "[CisplatinSTAD] rf_rep: $BASE/2rf_rep.R"
"$RBIN" "$BASE/2rf_rep.R"

echo "[CisplatinSTAD] done"
