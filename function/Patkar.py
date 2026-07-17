import os
import time
import numpy as np
import pandas as pd


def PatkarCal(exprSetFile, link_file, save_path, threshold=1.0):
    """
    Python reimplementation of Patkar et al. 2017.

    Reference:
        Patkar S et al., PLoS Comput Biol 13(11): e1005793, 2017.
        https://doi.org/10.1371/journal.pcbi.1005793

    Algorithm (faithful to original paper):
        For each sample, obtain a sample-specific PIN by removing nodes
        with expression < threshold. Edge weight = original STRING confidence
        score, unchanged across samples.

    Parameters
    ----------
    exprSetFile : str
        Expression matrix CSV (genes x samples, index_col=0).
    link_file : str
        PPI link CSV with columns [protein1, protein2, (score)].
    save_path : str
        Output directory. Each sample produces one file: {sample}.txt
        with columns [gene1, gene2, score].
    threshold : float
        Fixed expression cutoff. Genes below this are removed.
        Default = 1.0 (RPKM < 1 in the original paper).
    """

    os.makedirs(save_path, exist_ok=True)

    # ---------- Load expression matrix ----------
    print("[INFO] Loading expression matrix...")
    expr_df = pd.read_csv(exprSetFile, index_col=0)
    expr_df.index = expr_df.index.astype(str)
    expr_df.columns = expr_df.columns.astype(str)

    samples = expr_df.columns.tolist()
    gene_set = set(expr_df.index.tolist())
    print(f"[INFO] Samples: {len(samples)}, Genes: {len(gene_set)}")
    print(f"[INFO] Threshold: {threshold}")

    # ---------- Load PPI links ----------
    print("[INFO] Loading PPI links...")
    link_df = pd.read_csv(link_file, index_col=0)
    link_df['protein1'] = link_df['protein1'].astype(str).str.strip()
    link_df['protein2'] = link_df['protein2'].astype(str).str.strip()

    mask = (link_df['protein1'].isin(gene_set) &
            link_df['protein2'].isin(gene_set))
    link_df = link_df[mask].reset_index(drop=True)
    print(f"[INFO] PPI links after filtering: {len(link_df)}")

    if 'combined_score' in link_df.columns:
        w_col = 'combined_score'
    elif 'score' in link_df.columns:
        w_col = 'score'
    else:
        link_df['weight'] = 1.0
        w_col = 'weight'

    # ---------- Per-sample network construction ----------
    t_start = time.time()

    for s_idx, sample in enumerate(samples):
        print(f"[INFO] Sample {s_idx+1}/{len(samples)}: {sample}")

        expr_vals = expr_df[sample]
        expressed = set(expr_vals[expr_vals >= threshold].index.tolist())

        sample_links = link_df[
            link_df['protein1'].isin(expressed) &
            link_df['protein2'].isin(expressed)
        ].copy()

        if len(sample_links) == 0:
            print(f"[WARN]   No edges remain; skipping {sample}")
            continue

        df_out = sample_links[['protein1', 'protein2', w_col]].copy()
        df_out.columns = ['gene1', 'gene2', 'score']
        df_out.to_csv(
            os.path.join(save_path, f"{sample}.txt"),
            sep='\t', index=False
        )

    print(f"[INFO] Done. Elapsed: {time.time() - t_start:.1f}s")
