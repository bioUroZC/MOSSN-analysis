import os
import time
import numpy as np
import pandas as pd


def estimate_threshold(expr_matrix, method='percentile', percentile=25):
    """
    Estimate cohort-level expression threshold from the distribution of
    all non-zero values across the expression matrix.

    Parameters
    ----------
    expr_matrix : np.ndarray  (genes x samples)
    method : str  'percentile' | 'mean' | 'median'
    percentile : float  used when method='percentile'
    """
    flat = expr_matrix.flatten()
    nonzero = flat[flat > 0]
    if len(nonzero) == 0:
        return 0.0
    if method == 'percentile':
        return float(np.percentile(nonzero, percentile))
    elif method == 'mean':
        return float(np.mean(nonzero))
    elif method == 'median':
        return float(np.median(nonzero))
    else:
        raise ValueError(f"Unknown method: {method}")


def PPIXpressCal(exprSetFile, link_file, save_path,
                 threshold_method='percentile', percentile=25,
                 custom_threshold=None):
    """
    Python reimplementation of PPIXpress for gene-level expression data.

    Reference:
        Schaefer MH et al., Bioinformatics 2013.
        https://doi.org/10.1093/bioinformatics/btt441

    Note on simplification:
        The original PPIXpress operates at transcript level and uses
        domain-domain interaction data to assess isoform compatibility.
        Here we simplify to gene-level (consistent with how PPIXpress is
        used in ssPIN benchmarks), retaining its defining feature: a
        cohort-derived threshold rather than a fixed absolute cutoff.

    Algorithm (faithful to original design):
        1. Estimate a single threshold from the cohort expression matrix.
        2. For each sample, retain an edge only if BOTH genes are expressed
           above the threshold.
        3. Edge weight = original STRING confidence score, unchanged.

    Parameters
    ----------
    exprSetFile : str
        Expression matrix CSV (genes x samples, index_col=0).
    link_file : str
        PPI link CSV with columns [protein1, protein2, (score)].
    save_path : str
        Output directory. Each sample produces one file: {sample}.txt
        with columns [gene1, gene2, score].
    threshold_method : str
        'percentile', 'mean', or 'median'.
    percentile : float
        Used when threshold_method='percentile'. Default = 25.
    custom_threshold : float or None
        Override automatic estimation with a fixed value.
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

    # ---------- Estimate cohort threshold ----------
    if custom_threshold is not None:
        threshold = float(custom_threshold)
        print(f"[INFO] Using custom threshold: {threshold:.4f}")
    else:
        threshold = estimate_threshold(
            expr_df.values.astype(float),
            method=threshold_method,
            percentile=percentile
        )
        print(f"[INFO] Cohort threshold "
              f"({threshold_method}, p={percentile}): {threshold:.4f}")

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
