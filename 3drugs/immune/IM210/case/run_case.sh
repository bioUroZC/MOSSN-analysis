#!/bin/bash -l
#SBATCH -J case_IM210
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/3drugs/immune/logs/case_IM210_%j.out
#SBATCH -e /proj/c.zihao/work1/3drugs/immune/logs/case_IM210_%j.err

set -euo pipefail

DIR="/proj/c.zihao/work1/3drugs/immune/IM210/case"
mkdir -p /proj/c.zihao/work1/3drugs/immune/logs
cd "$DIR"

echo "[env] NODE=${SLURMD_NODENAME:-unknown} PART=${SLURM_JOB_PARTITION:-unknown} CPUS=${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

run_r() {
  echo "[run] $(date '+%F %T') $1"
  srun --cpu-bind=cores Rscript "$DIR/$1"
  echo "[done] $(date '+%F %T') $1"
}

run_r 01_define_immune_modules.R
run_r 02_compute_module_scores.R
run_r 03_response_trend.R
run_r 04_survival_analysis.R
run_r 05_visualize_case_sample.R
run_r 06_compare_ppi_gene.R
