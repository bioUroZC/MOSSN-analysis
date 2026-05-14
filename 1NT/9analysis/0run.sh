#!/bin/bash -l
#SBATCH -J atlas_pipeline
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/9analysis/logs/atlas_pipeline_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/9analysis/logs/atlas_pipeline_%j.err

set -euo pipefail

DIR="/proj/c.zihao/work1/1NT/9analysis"
mkdir -p "$DIR/logs" "$DIR/module_results" "$DIR/case_survival"
cd "$DIR"

echo "[env] NODE=${SLURMD_NODENAME:-unknown} PART=${SLURM_JOB_PARTITION:-unknown} CPUS=${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[run] $(date '+%F %T') $DIR/0runAtlasPipeline.R"
srun --cpu-bind=cores Rscript "$DIR/0runAtlasPipeline.R"
echo "[done] $(date '+%F %T') 0runAtlasPipeline.R"
