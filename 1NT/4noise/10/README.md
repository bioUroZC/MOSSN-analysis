# 3robust

This folder reruns the `2matrix` workflow on a noisy tumor-only `LUAD` expression matrix and compares the noisy results against the original `2matrix` full-sample results on the same tumor samples.

The default noise panel keeps:

- `MOSS_full`
- `MOSS_uniform`
- `SSN`
- `SWEET`
- `LIONESS`
- `Patkar`
- `PPIXpress`
- `Proteinarium`

The noise workflow keeps the external `Lung` normal reference unchanged and perturbs only the tumor expression matrix.

Suggested order:

1. `Rscript prepare.R`
2. `sbatch pipe.sh` or run `cal_*` scripts one by one
3. merge outputs with an R script analogous to `2matrix/3merge.R`
4. `python consistency.py`

Consistency definition:

- compare each noisy result against `2matrix/<Method>/merged_matrix.csv`
- use only tumor samples listed in `data/LUAD_noise_samples.csv`
- output:
  - `consistency/sample_level_consistency.csv`
  - `consistency/method_summary_consistency.csv`
  - `consistency/method_status.csv`

Batch run for all noise levels:

- `bash ../run_consistency.sh`
