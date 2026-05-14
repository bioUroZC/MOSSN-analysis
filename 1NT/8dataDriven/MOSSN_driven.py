import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix
def _rankdata_2d(X):
    temp = np.argsort(X, axis=1)
    ranks = np.empty_like(temp, dtype=float)
    idx = np.arange(X.shape[0])[:, None]
    ranks[idx, temp] = np.arange(1, X.shape[1] + 1).astype(float)
    return ranks


def prepare_data_MOSSN(expression_data):
    expression_data = expression_data.loc[expression_data.std(axis=1) > 0]

    # Spearman correlation via rank-transform then Pearson (memory-efficient)
    ranked = _rankdata_2d(expression_data.values).astype(float)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.sqrt((ranked ** 2).sum(axis=1, keepdims=True)) + 1e-10
    ranked /= norms
    corr_matrix = ranked @ ranked.T  # shape: (n_genes, n_genes)

    genes = expression_data.index.tolist()
    return corr_matrix, genes, expression_data


def build_graph_from_corr(corr_matrix, genes, expression_data, cor_threshold):
    rows_idx, cols_idx = np.where(
        np.triu(np.abs(corr_matrix) > cor_threshold, k=1)
    )

    G = nx.Graph()
    real_original_weights = {}
    for i, j in zip(rows_idx, cols_idx):
        u, v = genes[i], genes[j]
        w = float(abs(corr_matrix[i, j]))
        G.add_edge(u, v, weight=w)
        real_original_weights[(u, v)] = w

    used_genes = list(G.nodes)
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]

    return G, real_original_weights, expression_data


def MOSSN_single_sample(sample_id, G, real_original_weights, expression_data):
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):
            G[u][v]['weight'] = weight
            G[u][v]['adjusted_weight'] = weight
            G[u][v]['correction_score'] = 1.0

    sample_expression = expression_data[sample_id]
    sample_median = sample_expression.median()
    sample_iqr = sample_expression.quantile(0.75) - sample_expression.quantile(0.25)
    if sample_iqr == 0:
        sample_iqr = 1e-6

    expr_norm = ((sample_expression - sample_median) / sample_iqr).copy()
    alpha_mod = 1.0
    gamma = 2.0

    for u, v, data in G.edges(data=True):
        zu = expr_norm.get(u, 0)
        zv = expr_norm.get(v, 0)
        M_uv = zu + zv
        sigmoid = 1 / (1 + np.exp(-alpha_mod * M_uv))
        correction_score = 1 + gamma * (sigmoid - 0.5)
        correction_score = max(correction_score, 0)
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        data['adjusted_weight'] = BaseWeight * correction_score
        data['correction_score'] = correction_score

    nodes = list(G.nodes)
    n = len(nodes)
    node_to_index = {node: i for i, node in enumerate(nodes)}

    row, col, edge_weights = [], [], []
    for u, v in G.edges():
        weight = G[u][v]['adjusted_weight']
        row.append(node_to_index[u])
        col.append(node_to_index[v])
        edge_weights.append(weight)
        row.append(node_to_index[v])
        col.append(node_to_index[u])
        edge_weights.append(weight)

    T_sparse = csr_matrix((edge_weights, (row, col)), shape=(n, n))
    row_sums = np.array(T_sparse.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    T_sparse = T_sparse.multiply(1 / row_sums[:, None])

    q90 = sample_expression.quantile(0.9)
    seed_nodes = set(sample_expression[sample_expression > q90].index).intersection(G.nodes)

    rwr_alpha = 0.3
    P0 = np.array([1 if node in seed_nodes else 0 for node in nodes])
    if P0.sum() == 0:
        P0[:] = 1
    P0 = P0 / P0.sum()
    P = P0.copy()

    for _ in range(50):
        P_new = rwr_alpha * P0 + (1 - rwr_alpha) * T_sparse.T @ P
        if np.linalg.norm(P - P_new) < 1e-4:
            P = P_new
            break
        P = P_new

    ranks = np.argsort(np.argsort(P))
    P_normalized = ranks / (len(P) - 1)

    edges_data = []
    for u, v, data in G.edges(data=True):
        importance_u = P_normalized[node_to_index[u]]
        importance_v = P_normalized[node_to_index[v]]
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        final_weight = data['adjusted_weight'] * (importance_u + importance_v)

        edges_data.append({
            "Sample": sample_id,
            "Node1": u,
            "Node2": v,
            "BaseWeight": BaseWeight,
            "FinalWeight": final_weight,
        })

    edge_weights_df = pd.DataFrame(edges_data)
    edge_weights_df["FinalWeight"] = edge_weights_df["FinalWeight"].round(5)

    return edge_weights_df
