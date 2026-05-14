# 3robust

This folder reruns the `2matrix` workflow on a fixed subset of `LUAD` and compares the subset results against the original full-sample `2matrix` results.

The default robustness panel keeps:

- `MOSS_full`
- `MOSS_uniform`
- `SSN`
- `SWEET`
- `LIONESS`
- `Patkar`
- `PPIXpress`
- `Proteinarium`

Other MOSS ablation scripts are kept in this folder for separate ablation analyses, but they are not part of the default robustness comparison.

Suggested order:

1. `Rscript prepare.R`
2. `sbatch pipe.sh` or run `cal_*` scripts one by one
3. merge outputs with an R script analogous to `2matrix/3merge.R`
4. `python consistency.py`

Consistency definition:

- compare each subset result against `2matrix/<Method>/merged_matrix.csv`
- use only the selected tumor samples listed in `data/LUAD_half_samples.csv`
- output:
  - `consistency/sample_level_consistency.csv`
  - `consistency/method_summary_consistency.csv`
  - `consistency/method_status.csv`

Batch run for all ratios:

- `bash ../submit_robust_pipeline.sh`
