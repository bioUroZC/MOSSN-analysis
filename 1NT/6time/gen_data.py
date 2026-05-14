#!/usr/bin/env python3
"""Generate synthetic benchmark data for two timing dimensions.

Sample dimension  — fixed 10k edges, vary samples: 10, 20, 50, 100
  expr_S010.csv, expr_S020.csv, expr_S050.csv, expr_S100.csv
  links_E10000.csv

Network dimension — fixed 10 samples, vary edges: 10k, 20k, 50k, 100k
  expr_S010.csv   (reused from above)
  links_E10000.csv, links_E20000.csv, links_E50000.csv, links_E100000.csv

Normal reference uses the real GTEx file directly (Lung organ for SSN).
"""

import os
import numpy as np
import pandas as pd

SEED = 42
rng  = np.random.default_rng(SEED)

DATA_DIR = "/proj/c.zihao/work1/1NT/6time/data"
os.makedirs(DATA_DIR, exist_ok=True)

HURI_FILE    = "/proj/c.zihao/work1/1NT/1data/HuRI/links.csv"
SAMPLE_SIZES = [10, 20, 50, 100]
EDGE_COUNTS  = [10_000, 20_000, 50_000, 100_000]
N_MARKERS    = 200

# ── Gene pool ──────────────────────────────────────────────────────────────────
huri     = pd.read_csv(HURI_FILE, index_col=0)
all_genes = sorted(set(huri["protein1"]) | set(huri["protein2"]))
gene_arr  = np.array(all_genes)
n_genes   = len(gene_arr)
print(f"Gene pool: {n_genes} genes from HuRI")

# ── Expression matrices (one per sample size, all genes as rows) ──────────────
for N in SAMPLE_SIZES:
    n_tumor  = N // 2
    n_normal = N - n_tumor
    sample_ids = (
        [f"SYN_{i:03d}_01A" for i in range(1, n_tumor  + 1)] +
        [f"SYN_{i:03d}_11A" for i in range(1, n_normal + 1)]
    )
    expr = rng.lognormal(mean=2.0, sigma=1.0, size=(n_genes, N))
    pd.DataFrame(expr, index=all_genes, columns=sample_ids).to_csv(
        os.path.join(DATA_DIR, f"expr_S{N:03d}.csv"))
    print(f"  expr_S{N:03d}.csv  ({n_genes} genes × {N} samples)")

# ── Link files (one per edge count) ───────────────────────────────────────────
rng_links = np.random.default_rng(SEED)
max_edges = max(EDGE_COUNTS)

# Generate a large pool of unique random pairs, then subset
i_idx = rng_links.integers(0, n_genes, size=max_edges * 8)
j_idx = rng_links.integers(0, n_genes, size=max_edges * 8)
valid = i_idx != j_idx
i_idx, j_idx = i_idx[valid], j_idx[valid]
pairs = np.sort(np.stack([i_idx, j_idx], axis=1), axis=1)
_, first_occ = np.unique(pairs, axis=0, return_index=True)
pairs = pairs[np.sort(first_occ)]  # all unique pairs, in stable order

for K in EDGE_COUNTS:
    subset = pairs[:K]
    pd.DataFrame({
        "protein1": gene_arr[subset[:, 0]],
        "protein2": gene_arr[subset[:, 1]],
        "score":    1,
    }).to_csv(os.path.join(DATA_DIR, f"links_E{K:06d}.csv"))
    print(f"  links_E{K:06d}.csv  ({K} edges)")

# ── Markers (seed genes for MOSSN) ────────────────────────────────────────
marker_genes = rng.choice(gene_arr, size=N_MARKERS, replace=False)
pd.DataFrame(
    {"ensembl_gene_id": [f"ENSG{i:011d}" for i in range(N_MARKERS)],
     "symbol": marker_genes},
    index=[f"synthetic.{i+1}" for i in range(N_MARKERS)],
).to_csv(os.path.join(DATA_DIR, "markers.csv"))
print(f"  markers.csv  ({N_MARKERS} seed genes)")

print("gen_data.py done.")
