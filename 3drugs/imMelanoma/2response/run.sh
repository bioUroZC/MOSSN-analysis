#!/bin/bash -l
#SBATCH -J imMelanoma_response
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o /proj/c.zihao/work1/3drugs/logs/imMelanoma_response_%j.out
#SBATCH -e /proj/c.zihao/work1/3drugs/logs/imMelanoma_response_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/3drugs/imMelanoma/2response"

echo "[imMelanoma] prepare: $BASE/1prepare.R"
"$RBIN" "$BASE/1prepare.R"

echo "[imMelanoma] rf_dataset: $BASE/2rf_dataset.R"
"$RBIN" "$BASE/2rf_dataset.R"

echo "[imMelanoma] rf_rep: $BASE/2rf_rep.R"
"$RBIN" "$BASE/2rf_rep.R"

echo "[imMelanoma] done"
