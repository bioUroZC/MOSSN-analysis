#!/bin/bash -l
#SBATCH -J noise_submit
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/4coupled/noise/logs/noise_submit_%j.out
#SBATCH -e /proj/c.zihao/work1/4coupled/noise/logs/noise_submit_%j.err

set -euo pipefail

PYBIN="$HOME/.conda/envs/work1/bin/python"
RBIN="/usr/bin/Rscript"
NOISE_DIR="/proj/c.zihao/work1/4coupled/noise"
LOG_DIR="$NOISE_DIR/logs"
MODE="${MODE:-submit}"

mkdir -p "$LOG_DIR"

submit_jobs() {
    local prep_jobid
    local array_jobid
    local post_jobid

    prep_jobid=$(
        sbatch --parsable \
            --job-name=noise_prep \
            --nodes=1 \
            --cpus-per-task=8 \
            --mem=40G \
            -o "$LOG_DIR/noise_prep_%j.out" \
            -e "$LOG_DIR/noise_prep_%j.err" \
            --export=ALL,MODE=prep \
            "$0"
    )

    array_jobid=$(
        sbatch --parsable \
            --dependency=afterok:${prep_jobid} \
            --job-name=noise_test \
            --nodes=1 \
            --cpus-per-task=8 \
            --mem=40G \
            --array=0-5 \
            -o "$LOG_DIR/noise_%A_%a.out" \
            -e "$LOG_DIR/noise_%A_%a.err" \
            --export=ALL,MODE=compute \
            "$0"
    )

    post_jobid=$(
        sbatch --parsable \
            --dependency=afterok:${array_jobid} \
            --job-name=noise_post \
            --nodes=1 \
            --cpus-per-task=8 \
            --mem=40G \
            -o "$LOG_DIR/noise_post_%j.out" \
            -e "$LOG_DIR/noise_post_%j.err" \
            --export=ALL,MODE=post \
            "$0"
    )

    echo "[submit] noise prepare job: ${prep_jobid}"
    echo "[submit] compute array job: ${array_jobid}"
    echo "[submit] post-processing job: ${post_jobid} (afterok:${array_jobid})"
}

run_prep() {
    export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

    echo "[start] $(date '+%F %T')"
    echo "[env] NODE=$HOSTNAME CPUS=$SLURM_CPUS_PER_TASK"
    echo "[prep] generate noisy inputs"
    srun --cpu-bind=cores "$PYBIN" "$NOISE_DIR/0_add_noise.py"
    echo "[done] $(date '+%F %T')"
}

run_compute() {
    export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

    echo "[env] NODE=$HOSTNAME CPUS=$SLURM_CPUS_PER_TASK TASK=$SLURM_ARRAY_TASK_ID"

    SCRIPTS=(
        "$NOISE_DIR/function/cal_MOSSN_EXP_noise.py"
        "$NOISE_DIR/function/cal_MOSSN_NoCross_noise.py"
        "$NOISE_DIR/function/cal_MOSSN_Restart_noise.py"
        "$NOISE_DIR/function/cal_MOSSN_Direct_noise.py"
        "$NOISE_DIR/function/cal_MOSSN_DirectNoDyn_noise.py"
        "$NOISE_DIR/function/cal_MOSSN_MultiLayer_noise.py"
    )

    SCRIPT=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
    [ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

    echo "[compute] Running: $SCRIPT"
    srun --cpu-bind=cores "$PYBIN" "$SCRIPT"
}

run_post() {
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

    echo "[start] $(date '+%F %T')"
    echo "[env] NODE=$HOSTNAME CPUS=$SLURM_CPUS_PER_TASK"
    echo "[step] merge noise results"
    "$RBIN" "$NOISE_DIR/1_merge_noise.R"

    echo "[step] compute Spearman robustness"
    "$RBIN" "$NOISE_DIR/2_spearman_noise.R"

    echo "[step] plot robustness summary"
    "$RBIN" "$NOISE_DIR/3_plot_noise.R"

    echo "[done] $(date '+%F %T')"
}

case "$MODE" in
    submit)
        submit_jobs
        ;;
    prep)
        run_prep
        ;;
    compute)
        run_compute
        ;;
    post)
        run_post
        ;;
    *)
        echo "[ERR] Unknown MODE: $MODE" >&2
        exit 2
        ;;
esac
