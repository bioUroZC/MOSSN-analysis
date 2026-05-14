# 5parameter

Parameter sensitivity analysis for `LUAD` using `MOSSN_noPrior` with AUC-based evaluation on matched tumor-normal samples. The active parameter set is `gamma`, `restart`, and `seed_threshold`.

Modules:

- `prepare.R`: prepare paired `LUAD` samples (`01A` and `11A`)
- `gamma/`: sensitivity to `gamma`, with `rwr_alpha=0.3` and `seed_quantile=0.9`; AUC is summarized by `auc_gamma.R`
- `restart/`: sensitivity to restart probability for `MOSSN_noPrior`; AUC is summarized by `auc_restart.R`
- `seed_threshold/`: sensitivity to the high-expression seed quantile, with `rwr_alpha=0.3`; AUC is summarized by `auc_seed_threshold.R`

Recommended order:

1. `Rscript prepare.R`
2. Run `restart/cal_restart_mossn.py`
3. Run `restart/auc_restart.R`
4. Run `seed_threshold/cal_seed_threshold_mossn.py`
5. Run `seed_threshold/auc_seed_threshold.R`
6. Run `gamma/cal_gamma_mossn.py`
7. Run `gamma/auc_gamma.R`

SLURM scripts:

- `slurm_prepare.sh`: generate paired LUAD expression and metadata
- `restart/slurm_restart.sh`: run restart sensitivity matrix generation, AUC summary, and plots
- `seed_threshold/slurm_seed_threshold.sh`: run seed-threshold sensitivity matrix generation, AUC summary, and plots
- `slurm_5parameter.sh`: run the active three-parameter pipeline in one job

Suggested submissions:

1. `sbatch slurm_prepare.sh`
2. `sbatch restart/slurm_restart.sh`
3. `sbatch seed_threshold/slurm_seed_threshold.sh`
4. `sbatch gamma/slurm_gamma.sh`

Or submit one full job:

- `sbatch slurm_5parameter.sh`
