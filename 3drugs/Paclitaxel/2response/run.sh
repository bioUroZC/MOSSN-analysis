#!/bin/bash -l
#SBATCH -J Paclitaxel_response
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o /proj/c.zihao/work1/3drugs/logs/Paclitaxel_response_%j.out
#SBATCH -e /proj/c.zihao/work1/3drugs/logs/Paclitaxel_response_%j.err

set -euo pipefail

RBIN="${RBIN:-Rscript}"
BASE="/proj/c.zihao/work1/3drugs/Paclitaxel/2response"

echo "[Paclitaxel] prepare: $BASE/1prepare.R"
"$RBIN" "$BASE/1prepare.R"

echo "[Paclitaxel] rf_dataset: $BASE/2rf_dataset.R"
"$RBIN" "$BASE/2rf_dataset.R"

echo "[Paclitaxel] rf_rep: $BASE/2rf_rep.R"
"$RBIN" "$BASE/2rf_rep.R"

echo "[Paclitaxel] done"
