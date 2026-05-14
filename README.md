# mossn-analysis

Code for the MOSSN analysis project and its downstream evaluation workflows.

This repository is organized around four major analysis tracks:

- `function/`: core Python and R implementations of MOSSN and comparison methods
- `1NT/`: matched tumor-normal analyses, robustness studies, parameter sweeps, timing, and case studies
- `2survival/`: survival prediction workflows across multiple cancer cohorts
- `3drugs/`: drug response and immunotherapy response analyses
- `4coupled/`: coupled multi-omic extensions of MOSSN
- `0ref/`: reference-data preparation utilities

## Repository Notes

- This repository is prepared for code sharing before manuscript submission.
- Many workflows were originally run on an HPC cluster with `bash` or `sbatch`; the checked-in scripts preserve that execution style.
- Public upload should focus on code and lightweight documentation. Large intermediate outputs, downloaded datasets, and local runtime artifacts should stay out of Git history.

## Suggested Starting Points

- Core single-omic method: `function/MOSSN_noPrior.py`
- Coupled multi-omic method: `function/MOSSN_coupled.py`
- Parameter sensitivity example: `1NT/5parameter/README.md`
- Runtime benchmark example: `1NT/6time/README.md`
- Pan-cancer atlas writing outline: `1NT/9analysis/RESULTS_atlas_outline.md`

## Environment

The codebase mixes:

- Python
- R
- shell scripts for batch execution

Common Python dependencies used in the core methods include:

- `numpy`
- `pandas`
- `networkx`
- `scipy`

R package requirements vary by subproject and are not yet centralized in a lockfile.

## Before Uploading To GitHub

- Review whether any downloaded cohort data or generated figures still need to be removed.
- Add a license file if you already know the release license for the manuscript code.
- Initialize Git locally if you have not done so yet: `git init`

