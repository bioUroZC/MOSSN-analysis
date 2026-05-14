# MOSSN-analysis

Code for the MOSSN project and its downstream analyses, including matched tumor-normal network analysis, survival modeling, drug response analysis, and coupled multi-omic extensions.

## Overview

This repository is the analysis codebase used for the MOSSN study. It is shared primarily to document the computational workflow, analysis structure, and method implementation used in the manuscript.

The repository currently serves best as:

- a code archive for the project
- a reference for analysis logic and pipeline structure
- a starting point for readers who want to inspect or adapt the workflows

## Repository Structure

- `function/`: core Python and R implementations of MOSSN and comparison methods
- `0ref/`: reference-data preparation utilities
- `1NT/`: matched tumor-normal analyses, robustness studies, parameter sweeps, timing benchmarks, and case studies
- `2survival/`: survival prediction workflows across multiple cancer cohorts
- `3drugs/`: drug response and immunotherapy response analyses
- `4coupled/`: coupled multi-omic extensions of MOSSN

## Suggested Entry Points

- Core single-omic method: `function/MOSSN_noPrior.py`
- Coupled multi-omic method: `function/MOSSN_coupled.py`
- Parameter sensitivity example: `1NT/5parameter/README.md`
- Runtime benchmark example: `1NT/6time/README.md`
- Pan-cancer atlas writing outline: `1NT/9analysis/RESULTS_atlas_outline.md`

## Environment Notes

This codebase mixes:

- Python
- R
- shell scripts for batch execution

Common Python dependencies used in the core methods include:

- `numpy`
- `pandas`
- `networkx`
- `scipy`

R package requirements vary by subproject and are not yet centralized in a lockfile.

## Important Note On Paths

Many scripts in this repository were originally developed and executed in a personal HPC environment. For that reason, a substantial number of scripts still contain absolute paths such as `/proj/c.zihao/work1/...`.

For the current GitHub release, these paths are intentionally retained in many places because the repository is being shared mainly for code display and workflow reference rather than as a fully cleaned, turnkey reproduction package.

In practice, this means:

- the repository is useful for understanding the analysis workflow and method implementation
- some scripts will require local path adjustment before they can be executed in another environment
- not all workflows are expected to run out of the box on a new machine without adaptation

## Scope Of This Release

- Historical `old` analysis folders have been removed from the shared repository snapshot.
- Core dataset lists used by major survival, drug-response, and coupled-analysis pipelines were synchronized to match their active analysis scripts.
- The repository currently focuses on code and lightweight documentation rather than bundled result files or downloaded raw datasets.

## Notes For Readers

- If you are reading this repository alongside the manuscript, start from the folder most relevant to the corresponding analysis section.
- If you want to rerun a workflow, first inspect its local `README.md`, `1prepare.R`, `2merge.R`, or `run*.sh` entry scripts.
- If you want to adapt the code to a new machine, absolute paths are the first thing to replace.
