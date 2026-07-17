import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def MOSSN_noCorr_single_sample(
    sample_id,
    G,
    uniform_original_weights,
    expression_data,
    rwr_alpha=0.3,
    seed_quantile=0.9,
):
    # Ablation: no expression-based correction -> edges keep their base weight.
    for (u, v), weight in uniform_original_weights.items():
        if G.has_edge(u, v):
            G[u][v]["weight"] = weight
            G[u][v]["adjusted_weight"] = weight

    nodes = list(G.nodes)
    n = len(nodes)
    node_to_index = {node: i for i, node in enumerate(nodes)}

    row, col, edge_weights = [], [], []
    for u, v in G.edges():
        weight = G[u][v]["adjusted_weight"]
        row.append(node_to_index[u])
        col.append(node_to_index[v])
        edge_weights.append(weight)
        row.append(node_to_index[v])
        col.append(node_to_index[u])
        edge_weights.append(weight)

    t_sparse = csr_matrix((edge_weights, (row, col)), shape=(n, n))
    row_sums = np.array(t_sparse.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    t_sparse = t_sparse.multiply(1 / row_sums[:, None])

    sample_expression = expression_data[sample_id]
    seed_nodes = set(
        sample_expression[sample_expression > sample_expression.quantile(seed_quantile)].index
    ).intersection(G.nodes)
    p0 = np.array([1 if node in seed_nodes else 0 for node in nodes])
    if p0.sum() == 0:
        p0[:] = 1
    p0 = p0 / p0.sum()
    p = p0.copy()

    for _ in range(50):
        p_new = rwr_alpha * p0 + (1 - rwr_alpha) * t_sparse.T @ p
        if np.linalg.norm(p - p_new) < 1e-4:
            p = p_new
            break
        p = p_new

    # Average ranks so that tied probabilities (e.g. zero-probability nodes in
    # seedless components) receive identical importance instead of arbitrary
    # distinct ranks determined by node order.
    ranks = rankdata(p, method="average")
    p_normalized = (ranks - 1) / max(len(p) - 1, 1)

    edges_data = []
    for u, v in G.edges():
        importance_u = p_normalized[node_to_index[u]]
        importance_v = p_normalized[node_to_index[v]]
        base_weight = uniform_original_weights.get((u, v), uniform_original_weights.get((v, u), 0))
        final_weight = base_weight * (importance_u + importance_v)
        edges_data.append(
            {
                "Sample": sample_id,
                "Node1": u,
                "Node2": v,
                "BaseWeight": base_weight,
                "FinalWeight": final_weight,
            }
        )

    edge_weights_df = pd.DataFrame(edges_data)
    edge_weights_df["FinalWeight"] = edge_weights_df["FinalWeight"].round(5)
    return edge_weights_df
