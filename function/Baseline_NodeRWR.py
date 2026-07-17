import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def prepare_data_NodeRWR(links, expression_data):
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])
    links = link_filtered
    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]

    G = nx.Graph()
    for _, row in links.iterrows():
        G.add_edge(row['protein1'], row['protein2'], weight=row['score'])

    real_original_weights = {
        (row['protein1'], row['protein2']): row['score']
        for _, row in links.iterrows()
    }

    return G, real_original_weights, expression_data


def NodeRWR_single_sample(sample_id, G, real_original_weights, expression_data,
                           gamma=2.0, rwr_alpha=0.3, seed_quantile=0.9):
    """Same expression-correction + RWR steps as MOSSN_noPrior_single_sample,
    but returns per-node importance scores instead of aggregating them into
    edge scores. Used to test whether MOSSN's edge-level representation
    adds anything over the node-level RWR scores it is derived from.
    """
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

    seed_nodes = set(sample_expression[sample_expression > sample_expression.quantile(seed_quantile)].index).intersection(G.nodes)

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

    # Average ranks so that tied probabilities (e.g. zero-probability nodes in
    # seedless components) receive identical importance instead of arbitrary
    # distinct ranks determined by node order.
    ranks = rankdata(P, method="average")
    P_normalized = (ranks - 1) / max(len(P) - 1, 1)

    return pd.Series(P_normalized, index=nodes, name=sample_id)
