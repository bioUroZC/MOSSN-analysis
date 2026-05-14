#!/bin/bash
set -euo pipefail

BASE="/proj/c.zihao/work1/1NT/4noise"

mkdir -p "$BASE/logs"

jid_prepare=$(sbatch "$BASE/slurm_prepare_noise.sh" | awk '{print $4}')
jid_pipe=$(sbatch --dependency=afterok:"$jid_prepare" "$BASE/slurm_pipe_noise.sh" | awk '{print $4}')
jid_merge=$(sbatch --dependency=afterok:"$jid_pipe" "$BASE/slurm_merge_consistency.sh" | awk '{print $4}')
jid_plot=$(sbatch --dependency=afterok:"$jid_merge" "$BASE/slurm_plot_consistency.sh" | awk '{print $4}')

echo "prepare $jid_prepare"
echo "pipe $jid_pipe"
echo "merge_consistency $jid_merge"
echo "plot $jid_plot"
