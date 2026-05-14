import os
import time
import numpy as np
import pandas as pd
import networkx as nx


def _resolve_edge_costs(link_df, max_path_length, max_path_cost):
    """
    Infer how edge scores should be converted into Dijkstra costs.

    Supported cases:
        1. Weighted confidence scores above 1:
           cost = max(score) - score
        2. Normalized confidence scores in [0, 1]:
           cost = 1 - score
        3. Unweighted networks (all scores identical or no score column):
           cost = 1 for every edge
           and path cost falls back to path length control
    """

    if 'combined_score' in link_df.columns:
        score_col = 'combined_score'
    elif 'score' in link_df.columns:
        score_col = 'score'
    else:
        print("[INFO] No score column found; treating interactome as unweighted.")
        return np.ones(len(link_df), dtype=float), max_path_length

    score_series = pd.to_numeric(link_df[score_col], errors='coerce')
    valid_scores = score_series.dropna()

    if valid_scores.empty:
        print("[INFO] Score column is empty/non-numeric; treating interactome as unweighted.")
        return np.ones(len(link_df), dtype=float), max_path_length

    score_min = float(valid_scores.min())
    score_max = float(valid_scores.max())
    score_unique = int(valid_scores.nunique())

    if score_unique <= 1:
        print("[INFO] Scores are constant; treating interactome as unweighted.")
        return np.ones(len(link_df), dtype=float), max_path_length

    if 0.0 <= score_min and score_max <= 1.0:
        print("[INFO] Detected normalized edge scores in [0, 1]; using cost = 1 - score.")
        return 1.0 - score_series.to_numpy(dtype=float), max_path_cost

    print(f"[INFO] Detected weighted edge scores in [{score_min:.4g}, {score_max:.4g}]; "
          "using cost = max_score - score.")
    effective_max_path_cost = max_path_cost
    if max_path_cost <= max_path_length + 1e-12:
        effective_max_path_cost = max_path_cost * score_max
    return score_max - score_series.to_numpy(dtype=float), effective_max_path_cost


def ProteinariumCal(exprSetFile, link_file, save_path,
                    n_seed=30,
                    max_path_length=2,
                    max_path_cost=2.0):
    """
    Python reimplementation of Proteinarium for expression matrix input.

    Reference:
        Armanious D et al., Genomics 112(6): 4288-4296, 2020.
        https://doi.org/10.1016/j.ygeno.2020.07.028
        APPIC application note using Proteinarium-style construction:
        https://academic.oup.com/narcancer/article/7/1/zcae047/7954504

    Algorithm (faithful to original paper):
        1. For each sample, select seed genes (top n_seed by expression).
        2. Map seed genes onto STRING interactome.
        3. For each pair of seed proteins, run Dijkstra's shortest path.
           Edge cost is defined so that higher-confidence edges have
           lower cost.
        4. Retain only paths where:
               number of edges <= max_path_length  AND
               sum of edge costs  <= max_path_cost
        5. The union of all retained path edges forms the sample network.

    Output:
        Binary edge presence/absence per sample (1 = edge present, 0 = absent).
        Only edges present in the sample network are written to file (score=1).
        Absent edges are implicitly 0.

    Parameters
    ----------
    exprSetFile : str
        Expression matrix CSV (genes x samples, index_col=0).
    link_file : str
        PPI link CSV with columns [protein1, protein2, combined_score/score].
    save_path : str
        Output directory. Each sample produces {sample}.txt with columns
        [gene1, gene2, score] where score is always 1.
    n_seed : int
        Number of top-expressed genes used as seed genes per sample.
        Default = 30. Proteinarium used 10 seed genes in simulation and
        larger values such as 50-100 in expression-based use cases:
        https://www.sciencedirect.com/science/article/pii/S0888754320303050
    max_path_length : int
        Maximum number of edges in a retained path.
        Default = 2. This follows a later Proteinarium-based application
        (APPIC) that used maxPathLength = 2, i.e. allowing up to one
        intermediary node while limiting over-expansion of the network:
        https://academic.oup.com/narcancer/article/7/1/zcae047/7954504
    max_path_cost : float
        Maximum cumulative edge cost for a retained path.
        Default = 2.0. This is appropriate for normalized [0, 1] edge
        scores when using cost = 1 - score, and also matches the
        unweighted case when max_path_length = 2. For weighted scores
        above 1, this value is automatically rescaled when it is still
        given in normalized units.
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

    # ---------- Load PPI and build global graph ----------
    print("[INFO] Loading PPI links and building global graph...")
    link_df = pd.read_csv(link_file, index_col=0)
    link_df['protein1'] = link_df['protein1'].astype(str).str.strip()
    link_df['protein2'] = link_df['protein2'].astype(str).str.strip()

    mask = (link_df['protein1'].isin(gene_set) &
            link_df['protein2'].isin(gene_set))
    link_df = link_df[mask].reset_index(drop=True)
    print(f"[INFO] PPI links after filtering: {len(link_df)}")

    edge_costs, effective_max_path_cost = _resolve_edge_costs(
        link_df=link_df,
        max_path_length=max_path_length,
        max_path_cost=max_path_cost
    )
    print(f"[INFO] Effective max path cost: {effective_max_path_cost}")

    G_global = nx.Graph()
    for (_, row), cost in zip(link_df.iterrows(), edge_costs):
        g1, g2 = row['protein1'], row['protein2']
        G_global.add_edge(g1, g2, cost=cost)

    print(f"[INFO] Global graph: {G_global.number_of_nodes()} nodes, "
          f"{G_global.number_of_edges()} edges")

    # ---------- Per-sample network construction ----------
    t_start = time.time()

    for s_idx, sample in enumerate(samples):
        print(f"[INFO] Sample {s_idx+1}/{len(samples)}: {sample}")

        expr_vals = expr_df[sample]

        # Select seed genes: top n_seed by expression
        candidates = expr_vals[expr_vals.index.isin(G_global.nodes)]
        if len(candidates) == 0:
            print(f"[WARN]   No seed genes in graph; skipping {sample}")
            continue

        seed_genes = (candidates
                      .nlargest(min(n_seed, len(candidates)))
                      .index.tolist())

        # Dijkstra between all seed pairs
        sample_edges = set()

        for i in range(len(seed_genes)):
            for j in range(i + 1, len(seed_genes)):
                src, dst = seed_genes[i], seed_genes[j]

                if src not in G_global or dst not in G_global:
                    continue

                try:
                    path = nx.dijkstra_path(
                        G_global, src, dst, weight='cost'
                    )
                except nx.NetworkXNoPath:
                    continue

                path_edge_count = len(path) - 1
                if path_edge_count > max_path_length:
                    continue

                path_cost = sum(
                    G_global[path[k]][path[k+1]]['cost']
                    for k in range(len(path) - 1)
                )
                if path_cost > effective_max_path_cost:
                    continue

                for k in range(len(path) - 1):
                    u, v = path[k], path[k + 1]
                    sample_edges.add((min(u, v), max(u, v)))

        if len(sample_edges) == 0:
            print(f"[WARN]   No valid paths found; skipping {sample}")
            continue

        # Output: binary (score=1 for all present edges)
        rows = [[g1, g2, 1] for (g1, g2) in sample_edges]
        df_out = pd.DataFrame(rows, columns=['gene1', 'gene2', 'score'])
        df_out.to_csv(
            os.path.join(save_path, f"{sample}.txt"),
            sep='\t', index=False
        )

    print(f"[INFO] Done. Elapsed: {time.time() - t_start:.1f}s")
