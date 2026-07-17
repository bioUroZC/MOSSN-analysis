import numpy as np
import pandas as pd
import networkx as nx


def prepare_data_RandomBackbone(links, expression_data, seed=1):
    """Erdos-Renyi random PPI backbone: same node set and same number of
    edges as the STRING backbone, wired uniformly at random with no degree
    constraint. Edge weights are the real STRING score values, shuffled
    onto the random edge set, so the prior-weight scale matches the real
    analysis and only the topology differs.

    Returns (G, real_original_weights, expression_data), same shape as
    MOSSN_noPrior.prepare_data_MOSSN_noPrior, so MOSSN_noPrior_single_sample
    can be reused unchanged for the per-sample computation.
    """
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])
    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]

    genes = list(used_genes)
    n_edges = len(link_filtered)
    scores = link_filtered['score'].to_numpy()

    rng = np.random.default_rng(seed)
    random_graph = nx.gnm_random_graph(len(genes), n_edges, seed=seed)
    random_graph = nx.relabel_nodes(random_graph, {i: g for i, g in enumerate(genes)})

    shuffled_scores = rng.permutation(scores)
    real_original_weights = {}
    for (u, v), w in zip(random_graph.edges(), shuffled_scores[: random_graph.number_of_edges()]):
        random_graph[u][v]['weight'] = float(w)
        real_original_weights[(u, v)] = float(w)

    return random_graph, real_original_weights, expression_data
