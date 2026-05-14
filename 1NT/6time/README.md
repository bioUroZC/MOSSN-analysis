# 6time

This directory benchmarks the computational efficiency of `MOSSN_noPrior` against
comparison methods on two dimensions:

- sample size scaling: fix `n_edges = 10,000`, vary `n_samples = 10, 20, 50, 100`
- network size scaling: fix `n_samples = 10`, vary `n_edges = 10,000, 20,000, 50,000, 100,000`

The current pipeline measures end-to-end runtime for each method on the same
synthetic benchmark grids, then aggregates repeated runs into summary tables and
plots.

## Files

- `gen_data.py`: generate synthetic expression matrices, edge lists, and marker genes
- `run_timing.py`: run Python methods for one `(n_samples, n_edges)` grid point
- `run_timing_LIONESS.R`: run LIONESS for the same benchmark grid
- `aggregate_results.py`: merge per-grid CSV outputs into plot-ready tables
- `plot_timing.R`: summarise results and generate figures
- `submit_timing.sh`: submit the benchmark jobs on SLURM

## How To Run

1. Generate data and submit timing jobs:

```bash
bash submit_timing.sh
```

2. After the jobs finish, merge per-grid outputs:

```bash
python aggregate_results.py
```

3. Build summary tables and plots:

```bash
Rscript plot_timing.R
```

## Outputs

Generated files are written under this directory:

- `data/`: synthetic benchmark inputs
- `results/timing_S*_E*.csv`: per-grid Python timing results
- `results/timing_LIONESS_S*_E*.csv`: per-grid LIONESS timing results
- `results/timing_results_py.csv`: merged Python timing results
- `results/timing_results_LIONESS.csv`: merged LIONESS timing results
- `results/timing_summary_samples.csv`: summary table for fixed-edge sample scaling
- `results/timing_summary_edges.csv`: summary table for fixed-sample edge scaling
- `plots/timing_line_samples.pdf`: time per sample vs sample size
- `plots/timing_line_edges.pdf`: time per sample vs edge count
- `plots/timing_bar_N100.pdf`: time per sample comparison at `N = 100, E = 10k`
- `plots/timing_memory_N100.pdf`: Python peak memory comparison at `N = 100, E = 10k`

## Notes

- `MOSSN_noPrior` is timed with file loading included, matching the other methods'
  end-to-end timing scope.
- The memory column is recorded for Python methods via `tracemalloc`; LIONESS
  keeps `peak_memory_mb = NA`.
- `plot_timing.R` keeps the sample-scaling and edge-scaling summaries separate
  to avoid mixing different benchmark dimensions.
