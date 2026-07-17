# MOSSN-analysis

Analysis workflows accompanying the MOSSN study, including matched tumor–normal network analysis, survival prediction, immunotherapy-response analysis, pan-cancer analysis, and coupled multi-omic extensions.

## Overview

MOSSN constructs sample-specific weighted interaction networks by integrating molecular profiles with a reference protein–protein interaction (PPI) network.

This repository contains the analysis code used to generate the results reported in the manuscript, including data preprocessing, network construction, benchmarking, robustness analyses, survival modeling, immunotherapy-response analyses, and multi-omic extensions.

The standalone Python implementation of MOSSN is available at:

- https://github.com/bioUroZC/MOSSN

## Quick start

```bash
git clone https://github.com/bioUroZC/MOSSN-analysis.git
cd MOSSN-analysis

# 1. Point every script at your working root (see "Path configuration")
export MOSSN_ROOT=$(pwd)

# 2. Install dependencies
pip install -r requirements.txt
Rscript install_R_packages.R

# 3. Download and prepare the input data, then run a module in numbered order, e.g.
Rscript 1NT/1data/string/1linkPrepare.R      # build the STRING reference network
python  1NT/2string/benchmark/cal_MOSSN_uniform.py
Rscript 1NT/2string/benchmark/3merge.R
Rscript 1NT/2string/benchmark/4AUCcal.R
```

## MOSSN implementation

The primary single-omic implementation is `function/MOSSN_uniform.py`. For each sample it:

1. Takes the reference PPI network as the set of candidate interactions.
2. Initializes all retained reference edges with a base weight of 1.0.
3. Adjusts edge weights using sample-specific molecular deviations.
4. Propagates node information with a random walk with restart.
5. Combines the resulting node importance scores with the adjusted edge weights to produce the final sample-specific network.

The coupled multi-omic implementation is `function/MOSSN_coupled.py`.

### Default parameters

The default parameters are identical across every analysis in this repository:

| Parameter | Value | Meaning |
| --- | --- | --- |
| `gamma` | 2.0 | Edge-correction strength |
| `alpha_mod` | 1.0 | Sigmoid sensitivity |
| `rwr_alpha` | 0.3 | RWR restart probability |
| `seed_quantile` | 0.9 | Seed-node threshold (90th percentile) |
| `tol` | 1e-4 | RWR convergence threshold |
| `max_iter` | 50 | RWR iteration cap |

## Repository structure

- `function/`: implementations of MOSSN, coupled MOSSN, comparison methods, and ablation variants
- `0ref/`: GTEx reference-expression preparation
- `1NT/`: matched tumor–normal analysis, benchmarking, robustness evaluation, parameter sensitivity, runtime analysis, and pan-cancer analyses. Reference-network preparation lives under `1NT/1data/` (`string/`, `HuRI/`, `biogrid/`, `intact/`)
- `2survival/`: leave-one-dataset-out survival prediction analyses
- `3drugs/`: immunotherapy-response analyses (IMvigor210 anti-PD-L1 cohort)
- `4coupled/`: coupled multi-omic extensions of MOSSN

## Suggested entry points

- Primary single-omic MOSSN implementation: `function/MOSSN_uniform.py`
- Coupled multi-omic implementation: `function/MOSSN_coupled.py`
- Reference-network preparation: `1NT/1data/string/1linkPrepare.R`
- Matched tumor–normal benchmarking: `1NT/2string/benchmark/`
- Ablation analyses: `1NT/2string/ablation/`
- Pan-cancer analyses: `1NT/9analysis/`
- Survival prediction workflows: `2survival/`
- Immunotherapy case study: `3drugs/immune/IM210/case/`

## Survival modeling

Survival prediction is evaluated using a leave-one-dataset-out strategy. For each held-out dataset:

1. Feature filtering and feature selection are performed using the training datasets only.
2. Model parameters are estimated using the training datasets.
3. Model performance is evaluated exclusively in the held-out dataset.
4. C-index and time-dependent AUC are reported for each held-out dataset.

Cancer-specific implementations are available under:

- `2survival/<cancer_type>/2quick/ml.R`

## Path configuration

All scripts resolve file locations from a single root directory, exposed as the environment variable `MOSSN_ROOT`. Every script defines it near the top:

```python
# Python
PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")
```

```r
# R
PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")
```

To run the workflows on another system, set `MOSSN_ROOT` once. No per-script editing is required:

```bash
export MOSSN_ROOT=/path/to/mossn-analysis
```

or, from within an R session:

```r
Sys.setenv(MOSSN_ROOT = "/path/to/mossn-analysis")
```

The fallback value `/proj/c.zihao/work1` is the HPC location where the reported analyses were originally executed. It is retained as a default so that the published code is exactly the code that produced the results; setting `MOSSN_ROOT` overrides it everywhere.

The repository layout mirrors the expected working-directory layout: scripts read and write under `$MOSSN_ROOT/1NT/`, `$MOSSN_ROOT/2survival/`, and so on. Downloaded and intermediate data are therefore placed alongside the corresponding scripts.

Changing `MOSSN_ROOT` does not alter the computational method or analysis parameters.

## Software requirements

The analysis workflows use Python and R.

The results reported in the manuscript were produced with **R 4.3.1** and **Python 3.10**.

### Python

Install with `pip install -r requirements.txt`. The dependencies are:

- `numpy`
- `pandas`
- `networkx`
- `scipy`
- `tqdm`

### R

Install with `Rscript install_R_packages.R`, which installs the required CRAN and Bioconductor packages. The required packages are also declared near the beginning of the corresponding R scripts.

## Input data and reference networks

The repository does not redistribute controlled-access or third-party datasets. Data sources and accession identifiers are described in the manuscript and supplementary materials.

Users should download the corresponding datasets and prepare them using the download and preprocessing scripts within each analysis module, for example:

- `1NT/1data/TCGA/1download.R` — TCGA expression and clinical data
- `1NT/1data/string/1linkPrepare.R` — STRING reference network
- `0ref/GTEx/1Prepare.R` — GTEx reference expression
- `3drugs/immune/IM210/1codes/1download.R` — IMvigor210 cohort preprocessing

The IMvigor210 anti-PD-L1 cohort is not redistributed here. Obtain the expression and response tables from the original publication's data package (Mariathasan et al., *Nature* 2018) and place them in `3drugs/immune/IM210/data/` as:

- `IMvigor210_exprSet.csv` — genes × samples expression matrix
- `IMvigor210_FollowUp.csv` — per-sample response and survival annotation

`1codes/1download.R` then reads these two files and writes the harmonised `exprSet.csv` and `pd.csv` used by the rest of the module.

The reference PPI network should contain at least the following columns:

- `protein1`
- `protein2`

All retained interactions are assigned an initial edge weight of 1.0.

## Reproducing individual analyses

To reproduce an analysis:

1. Set `MOSSN_ROOT` (see "Path configuration").
2. Identify the corresponding analysis directory.
3. Run the download and preprocessing scripts to obtain the required input data.
4. Run the scripts in their numbered order when numbered workflow files are provided. A typical module proceeds as: download → prepare → `cal_*` (per-sample network construction) → `*merge*` (assemble the edge-by-sample matrix) → analysis and plotting.
5. Random seeds are fixed in the scripts; no further action is required.

Large intermediate network matrices are not stored in the repository because of their size. They can be regenerated from the processed molecular data and reference interaction network using the supplied scripts.

## Method variants

Alongside the primary MOSSN implementation, `function/` contains the benchmark comparison methods and several ablation and control variants (for example `MOSSN_noCorr`, `EdgeNoRWR`, `NodeRWR`). These are used only for the benchmarking and ablation analyses.

## Reproducibility and versioning

The results in the manuscript correspond to the code version associated with the manuscript release. Users are encouraged to use the corresponding GitHub release or commit rather than the continuously updated development branch.

Random seeds are fixed wherever stochastic procedures are involved:

- R workflows call `set.seed()` before random sampling, k-means clustering, cross-validation, and network randomization. The seeds used are `1`, `42`, `123`, and `1234`, depending on the module.
- Python workflows derive seeds deterministically rather than drawing them at run time: `cal_6_RandomBackbone.py` seeds each dataset by its index, and `cal_8_PermutedControl.py` derives a per-sample seed from the sample identifier via `Baseline_PermutedControl.seed_from_sample_id(..., base_seed=1)`.

Given the same input data, the workflows are therefore deterministic.

## License

Released under the MIT License. See the `LICENSE` file for reuse conditions.

## Citation

If you use MOSSN or the accompanying analysis workflows, please cite the MOSSN manuscript and the software repository.
