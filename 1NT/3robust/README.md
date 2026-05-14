# 3robust

Robustness analysis for `LUAD` under four sample-retention ratios:

- `10`: keep 10% of tumor samples and 10% of Lung normal reference samples
- `30`: keep 30%
- `50`: keep 50%
- `70`: keep 70%

Each ratio follows the same workflow:

1. `prepare.R`: sample tumor and Lung normal subsets
2. `pipe.sh`: run 8 methods on the subset expression matrix
3. `merge.R`: convert raw outputs into `merged_matrix.csv`
4. `consistency.R`: compare subset results against the reference outputs from `../2matrix`
5. `plot_consistency.R`: summarize all ratios at the top level

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

- `slurm_prepare_robust.sh`: prepare subset inputs for `10/30/50/70`
- `slurm_pipe_robust.sh`: run all methods for all ratios as one array job
- `slurm_merge_consistency.sh`: run `merge.R` and `consistency.R` for each ratio
- `slurm_plot_consistency.sh`: generate top-level summary tables and PDF plots
- `submit_robust_pipeline.sh`: submit the full dependency chain automatically

Recommended submission:

`bash submit_robust_pipeline.sh`
