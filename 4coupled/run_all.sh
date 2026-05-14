#!/bin/bash
# Full pipeline: method generation -> 5merge -> 6stratification
# Usage: bash run_all.sh

set -euo pipefail

LOGDIR="/proj/c.zihao/work1/4coupled/logs"
mkdir -p "$LOGDIR"

# ── Step 1: compute the final 8 methods (array job) ───────────────────────
JID1=$(sbatch --parsable << 'EOF'
#!/bin/bash -l
#SBATCH -J 4coupled
#SBATCH -p GPU-A40
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=0-7
#SBATCH -o /proj/c.zihao/work1/4coupled/logs/4coupled_%A_%a.out
#SBATCH -e /proj/c.zihao/work1/4coupled/logs/4coupled_%A_%a.err

set -euo pipefail

PYBIN="$HOME/.conda/envs/work1/bin/python"
DIR="/proj/c.zihao/work1/4coupled/function"

export LD_LIBRARY_PATH="$HOME/.conda/envs/work1/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[env] NODE=$HOSTNAME CPUS=$SLURM_CPUS_PER_TASK"
"$PYBIN" -V

SCRIPTS=(
    "$DIR/cal_MOSSN_EXP.py"
    "$DIR/cal_MOSSN_MET.py"
    "$DIR/cal_MOSSN_CNV.py"
    "$DIR/cal_MOSSN_NoCross.py"
    "$DIR/cal_MOSSN_Restart.py"
    "$DIR/cal_MOSSN_Direct.py"
    "$DIR/cal_MOSSN_DirectNoDyn.py"
    "$DIR/cal_MOSSN_MultiLayer.py"
)

SCRIPT=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
[ -f "$SCRIPT" ] || { echo "[ERR] Script not found: $SCRIPT"; exit 2; }

echo "[task ${SLURM_ARRAY_TASK_ID}] Running: $SCRIPT"
srun --cpu-bind=cores "$PYBIN" "$SCRIPT"
EOF
)
echo "[1/3] Submitted array job: $JID1 (tasks 0-7)"

# ── Step 2: merge all method outputs into matrices ─────────────────────────
JID2=$(sbatch --parsable --dependency=afterok:$JID1 << 'EOF'
#!/bin/bash -l
#SBATCH -J 5merge
#SBATCH -p GPU-A40
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -o /proj/c.zihao/work1/4coupled/logs/merge_%j.out
#SBATCH -e /proj/c.zihao/work1/4coupled/logs/merge_%j.err

set -euo pipefail

RBIN="/usr/bin/Rscript"
echo "[env] NODE=$HOSTNAME MEM=40G"
"$RBIN" --version | head -1
srun "$RBIN" /proj/c.zihao/work1/4coupled/5merge.R
EOF
)
echo "[2/3] Submitted merge job: $JID2 (depends on $JID1)"

# ── Step 3: survival stratification + HR plots ─────────────────────────────
JID3=$(sbatch --parsable --dependency=afterok:$JID2 << 'EOF'
#!/bin/bash -l
#SBATCH -J surv_strat
#SBATCH -p GPU-A40
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH -o /proj/c.zihao/work1/4coupled/logs/surv_strat_%j.out
#SBATCH -e /proj/c.zihao/work1/4coupled/logs/surv_strat_%j.err

set -euo pipefail

RBIN="/usr/bin/Rscript"
echo "[env] NODE=$HOSTNAME CPUS=${SLURM_CPUS_PER_TASK:-NA} MEM=60G"
"$RBIN" --version | head -1

echo "[1/2] Running survival stratification"
srun "$RBIN" /proj/c.zihao/work1/4coupled/6stratification.R

echo "[2/2] Running HR plotting"
srun "$RBIN" /proj/c.zihao/work1/4coupled/7plot_survival_hr.R

echo "[done] Results under /proj/c.zihao/work1/4coupled/results/survival_stratification"
EOF
)
echo "[3/3] Submitted stratification job: $JID3 (depends on $JID2)"

echo ""
echo "Pipeline submitted. Check status with:"
echo "  squeue -j $JID1,$JID2,$JID3"
