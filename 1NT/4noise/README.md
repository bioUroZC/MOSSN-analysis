# 4noise

Noise robustness analysis for `LUAD` tumor samples under three additive noise levels:

- `05`: 5% gene-wise noise
- `10`: 10% gene-wise noise
- `20`: 20% gene-wise noise

Each level follows the same workflow:

1. `prepare.R`: add Gaussian noise to tumor samples only and keep the Lung normal reference unchanged
2. `pipe.sh`: run 8 methods on the noisy matrix
3. `merge.R`: convert raw method outputs into `merged_matrix.csv`
4. `consistency.R`: compare noisy results against the reference outputs from `../2matrix`
5. `plot_consistency.R`: summarize all noise levels at the top level

Methods:

- `MOSS_full`
- `MOSS_uniform`
- `SSN`
- `SWEET`
- `LIONESS`
- `Patkar`
- `PPIXpress`
- `Proteinarium`

Top-level SLURM entry points:

- `slurm_prepare_noise.sh`: prepare noisy input data for `05/10/20`
- `slurm_pipe_noise.sh`: run all methods for all noise levels as one array job
- `slurm_merge_consistency.sh`: run `merge.R` and `consistency.R` for each noise level
- `slurm_plot_consistency.sh`: generate top-level summary tables and PDF plots
- `submit_noise_pipeline.sh`: submit the full dependency chain automatically

Recommended submission:

`bash submit_noise_pipeline.sh`

This submits the jobs with dependencies:

1. `slurm_prepare_noise.sh`
2. `slurm_pipe_noise.sh` after prepare succeeds
3. `slurm_merge_consistency.sh` after pipe succeeds
4. `slurm_plot_consistency.sh` after merge/consistency succeeds

Manual submission is still possible if needed:

1. `sbatch slurm_prepare_noise.sh`
2. `sbatch --dependency=afterok:<prepare_jobid> slurm_pipe_noise.sh`
3. `sbatch --dependency=afterok:<pipe_jobid> slurm_merge_consistency.sh`
4. `sbatch --dependency=afterok:<merge_jobid> slurm_plot_consistency.sh`
