import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix


def prepare_data_MOSSN_noCorr(links, expression_data):
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


def MOSSN_noCorr_single_sample(sample_id, G, real_original_weights, expression_data,
                               rwr_alpha=0.3, seed_quantile=0.9):
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):
            G[u][v]['weight'] = weight
            G[u][v]['adjusted_weight'] = weight

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

    sample_expression = expression_data[sample_id]
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

    ranks = np.argsort(np.argsort(P))
    P_normalized = ranks / (len(P) - 1)

    edges_data = []
    for u, v in G.edges():
        importance_u = P_normalized[node_to_index[u]]
        importance_v = P_normalized[node_to_index[v]]
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        final_weight = BaseWeight * (importance_u + importance_v)

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
