# MOSSN-analysis

Analysis code accompanying:

> **MOSSN: A Multi-Omics Framework for Sample-Specific Protein Interaction Networks in Cancer**

MOSSN (Multi-Omics Single-Sample random walk for Network reweighting) assigns sample-specific weights to the edges of a reference protein-protein interaction network. It combines within-sample molecular deviations with random-walk propagation and supports expression, DNA methylation, and copy-number data.

This repository contains the scripts for benchmarking, ablation, robustness, survival, multi-omics, pan-cancer, and immunotherapy-response analyses.

## Method overview

For each sample, MOSSN:

1. Uses a reference PPI network as the candidate interaction backbone.
2. Modulates edge weights using within-sample molecular values.
3. Propagates molecular signals with a random walk with restart.
4. Combines corrected edge weights with propagated node importance scores.

The multi-omics extension uses the expression network as the primary layer and couples methylation and copy-number signals through gene-specific cross-layer edges.

## Repository structure

| Directory | Description |
| --- | --- |
| `function/` | Core MOSSN, multi-omics, benchmark, and ablation implementations. |
| `01data/` | Input-data and reference-network preparation. |
| `01NT/` | Tumor-normal, benchmark, ablation, topology, and GEO analyses. |
| `02analysis/` | Robustness, noise, network, runtime, and parameter analyses. |
| `03survival/` | Cancer-specific survival workflows. |
| `04multi/` | Multi-omics preparation, survival, and noise analyses. |
| `05atlas/` | Pan-cancer recurrent-edge, module, enrichment, and plotting analyses. |
| `06case/` | IMvigor210 immunotherapy-response case study. |

Most modules use numbered scripts in an approximate order such as:

```text
prepare/define -> calculate -> merge -> evaluate/plot
```

## Core implementations

- `function/ablation/MOSSN_func.py`: single-omics MOSSN functions.
- `function/multi/MOSSN_coupled.py`: coupled and multilayer multi-omics functions.
- `function/benchmark/`: SSN, SWEET, PPIXpress, Patkar, Proteinarium, and LIONESS implementations.
- `function/ablation/`: ablation and control models.

## Requirements

The workflows use Python and R.

Common Python packages include:

```text
numpy
pandas
networkx
scipy
tqdm
```

The R workflows use packages from the tidyverse, network-analysis, survival-analysis, and Bioconductor ecosystems.

## Getting started

1. Install Python, R, and the packages required by the selected workflow.
2. Prepare the input expression, clinical, survival, and reference-network files for that workflow.
3. Set the workflow-specific input and output paths.
4. Run the numbered scripts in the selected module in order.

## Running an analysis

Python scripts can be run with `python`; R scripts can be run with `Rscript`.

Representative entry points include:

```text
01data/06string/1linkPrepare.R
01NT/01string/02benchmark/01cal/
01NT/01string/01ablation/01cal/
02analysis/03net/
03survival/<cancer_type>/
04multi/
05atlas/01Pipeline.R
06case/02cal/
```

## Outputs

Depending on the module, the workflows generate sample-specific edge tables, edge-by-sample matrices, benchmark metrics, topology summaries, robustness results, survival features, multi-omics stratification results, pan-cancer summaries, and immunotherapy-response analyses.

## Reproducibility notes

- Reference networks are analyzed in separate STRING, HuRI, BioGRID, and IntAct workflows.
- MOSSN reweights a selected reference network; it does not infer a new interactome from scratch.
- Please consult the manuscript and supplementary materials for dataset definitions, preprocessing, parameters, and evaluation details.

## Citation

Please cite the accompanying MOSSN manuscript and the original sources of all datasets and reference networks used.
