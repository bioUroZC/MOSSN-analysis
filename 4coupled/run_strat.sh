#!/bin/bash
# Run only: 6stratification.R -> 7plot_survival_hr.R -> 8plot_natural_hr.R
# Usage: bash run_strat.sh

set -euo pipefail

LOGDIR="/proj/c.zihao/work1/4coupled/logs"
mkdir -p "$LOGDIR"

JID=$(sbatch --parsable << 'EOF'
#!/bin/bash -l
#SBATCH -J surv_strat
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH -o /proj/c.zihao/work1/4coupled/logs/surv_strat_%j.out
#SBATCH -e /proj/c.zihao/work1/4coupled/logs/surv_strat_%j.err

set -euo pipefail

RBIN="/usr/bin/Rscript"
echo "[env] NODE=$HOSTNAME CPUS=${SLURM_CPUS_PER_TASK:-NA} MEM=60G"
"$RBIN" --version | head -1

echo "[1/3] Running survival stratification"
srun "$RBIN" /proj/c.zihao/work1/4coupled/6stratification.R

echo "[2/3] Running balanced HR plotting"
srun "$RBIN" /proj/c.zihao/work1/4coupled/7plot_survival_hr.R

echo "[3/3] Running natural HR plotting"
srun "$RBIN" /proj/c.zihao/work1/4coupled/8plot_natural_hr.R

echo "[done] Results under /proj/c.zihao/work1/4coupled/results/survival_stratification"
EOF
)
echo "Submitted stratification job: $JID"

echo ""
echo "Check status with:"
echo "  squeue -j $JID"
