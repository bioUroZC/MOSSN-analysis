import numpy as np
import pandas as pd
import networkx as nx


def prepare_data_MOSSN_noRWR(links, expression_data):
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


def MOSSN_noRWR_single_sample(sample_id, G, real_original_weights, expression_data,
                              gamma=2.0):
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

    edges_data = []
    for u, v, data in G.edges(data=True):
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        edges_data.append({
            "Sample": sample_id,
            "Node1": u,
            "Node2": v,
            "BaseWeight": BaseWeight,
            "FinalWeight": data['adjusted_weight'],
        })

    edge_weights_df = pd.DataFrame(edges_data)
    edge_weights_df["FinalWeight"] = edge_weights_df["FinalWeight"].round(5)

    return edge_weights_df
