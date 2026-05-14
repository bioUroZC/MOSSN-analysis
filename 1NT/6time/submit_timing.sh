#!/bin/bash -l
# Step 1: generate data (run once, locally or as a short job)
# Step 2: sbatch this file to submit 7 Python timing jobs plus 1 R job

PYBIN="$HOME/.conda/envs/work1/bin/python"
DIR="/proj/c.zihao/work1/1NT/6time"
export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

mkdir -p "$DIR/logs" "$DIR/data" "$DIR/results" "$DIR/plots"

# ── Step 1: generate synthetic data ───────────────────────────────────────────
echo "Generating synthetic data..."
"$PYBIN" "$DIR/gen_data.py"

# ── Step 2: Python timing jobs (array of 7) ───────────────────────────────────
# Task → (N_SAMPLES, N_EDGES)
#   0: sample dim, N=10,  E=10k
#   1: sample dim, N=20,  E=10k
#   2: sample dim, N=50,  E=10k
#   3: sample dim, N=100, E=10k
#   4: network dim, N=10, E=20k
#   5: network dim, N=10, E=50k
#   6: network dim, N=10, E=100k
# LIONESS is handled separately below.

sbatch <<'PYEOF'
#!/bin/bash -l
#SBATCH -J timing_py
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --array=0-6
#SBATCH -o /proj/c.zihao/work1/1NT/6time/logs/timing_py_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/1NT/6time/logs/timing_py_%A_%a.err

PYBIN="$HOME/.conda/envs/work1/bin/python"
DIR="/proj/c.zihao/work1/1NT/6time"
export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

N_SAMPLES_ARR=(  10  20   50  100  10   10    10)
N_EDGES_ARR=(10000 10000 10000 10000 20000 50000 100000)

N=${N_SAMPLES_ARR[$SLURM_ARRAY_TASK_ID]}
K=${N_EDGES_ARR[$SLURM_ARRAY_TASK_ID]}

echo "[task ${SLURM_ARRAY_TASK_ID}] N_SAMPLES=$N  N_EDGES=$K"
srun "$PYBIN" "$DIR/run_timing.py" "$N" "$K"
PYEOF

# ── Step 3: LIONESS (R job, separate) ─────────────────────────────────────────
sbatch <<'REOF'
#!/bin/bash -l
#SBATCH -J timing_LIONESS
#SBATCH -N 1
#SBATCH -p GPU-A40
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/1NT/6time/logs/timing_LIONESS_%j.out
#SBATCH -e /proj/c.zihao/work1/1NT/6time/logs/timing_LIONESS_%j.err

DIR="/proj/c.zihao/work1/1NT/6time"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun Rscript "$DIR/run_timing_LIONESS.R"
REOF

echo "Jobs submitted."
echo "After the jobs finish, run:"
echo "  $PYBIN $DIR/aggregate_results.py"
echo "  Rscript $DIR/plot_timing.R"
