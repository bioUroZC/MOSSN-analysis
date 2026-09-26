"""Generate synthetic benchmark data for two timing dimensions.

Sample dimension  — fixed 10k edges, vary samples: 10, 20, 50, 100
  expr_S010.csv, expr_S020.csv, expr_S050.csv, expr_S100.csv
  links_E10000.csv

Network dimension — fixed 10 samples, vary edges: 10k, 20k, 50k, 100k
  expr_S010.csv   (reused from above)
  links_E10000.csv, links_E20000.csv, links_E50000.csv, links_E100000.csv

Gene names are taken from the STRING network; the edges themselves are random
pairs drawn over that gene pool, so edge count is controlled independently.
This stage creates only synthetic tumour/normal-labelled inputs.  The separate
paired-SSN timing runner supplies its cohort-native TCGA-LUAD normal reference.
"""

import os
import numpy as np
import pandas as pd

SEED = 42
rng  = np.random.default_rng(SEED)

DATA_DIR = "/proj/c.zihao/work1/02analysis/04time/data"
os.makedirs(DATA_DIR, exist_ok=True)

STRINGNET_FILE = "/proj/c.zihao/work1/01data/06string/links.csv"
SAMPLE_SIZES   = [10, 20, 50, 100]
EDGE_COUNTS    = [10_000, 20_000, 50_000, 100_000]

stringNet = pd.read_csv(STRINGNET_FILE, index_col=0)
all_genes = sorted(set(stringNet["protein1"]) | set(stringNet["protein2"]))
gene_arr  = np.array(all_genes)
n_genes   = len(gene_arr)
print(f"Gene pool: {n_genes} genes from STRING")

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

rng_links = np.random.default_rng(SEED)
max_edges = max(EDGE_COUNTS)

i_idx = rng_links.integers(0, n_genes, size=max_edges * 8)
j_idx = rng_links.integers(0, n_genes, size=max_edges * 8)
valid = i_idx != j_idx
i_idx, j_idx = i_idx[valid], j_idx[valid]
pairs = np.sort(np.stack([i_idx, j_idx], axis=1), axis=1)
_, first_occ = np.unique(pairs, axis=0, return_index=True)
pairs = pairs[np.sort(first_occ)]  

for K in EDGE_COUNTS:
    subset = pairs[:K]
    pd.DataFrame({
        "protein1": gene_arr[subset[:, 0]],
        "protein2": gene_arr[subset[:, 1]],
        "score":    1,
    }).to_csv(os.path.join(DATA_DIR, f"links_E{K:06d}.csv"))
    print(f"  links_E{K:06d}.csv  ({K} edges)")

print("01prepare.py done.")
