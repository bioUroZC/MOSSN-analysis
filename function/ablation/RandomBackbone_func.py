import pandas as pd
import networkx as nx

def RandomBackbone_prepare_data(links, expression_data, seed=1):
    """Erdos-Renyi random PPI backbone: same node set and same number of
    edges as the STRING backbone, wired uniformly at random with no degree
    constraint. Edge weights use the same uniform baseline as MOSSN, so only
    the topology differs from the real analysis.

    Returns (G, real_original_weights, expression_data), same shape as
    MOSSN_func.MOSSN_prepare_data, so MOSSN_single_sample can be reused
    unchanged for the per-sample computation.
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

    random_graph = nx.gnm_random_graph(len(genes), n_edges, seed=seed)
    random_graph = nx.relabel_nodes(random_graph, {i: g for i, g in enumerate(genes)})

    real_original_weights = {}
    for u, v in random_graph.edges():
        random_graph[u][v]['weight'] = 1.0
        real_original_weights[(u, v)] = 1.0

    return random_graph, real_original_weights, expression_data

def RandomBackbone_prepare_shared(links, seed=1):
    """Erdos-Renyi random PPI backbone built from the link table alone.

    RandomBackbone_prepare_data() draws its node set from one dataset's
    expression matrix, so two datasets with different gene coverage get
    different random graphs even under the same seed. Leave-one-dataset-out
    validation then has no edge that is present in every training dataset,
    and the whole feature set collapses. This variant takes the node set and
    the edge count from the link table, so one call produces a single
    topology that every dataset of a cancer can share.

    Returns (G, real_original_weights); the caller restricts its own
    expression matrix to G's nodes.
    """
    link_filtered = links.drop_duplicates(subset=['protein1', 'protein2'])
    genes = list(pd.unique(link_filtered[['protein1', 'protein2']].values.ravel()))
    n_edges = len(link_filtered)

    random_graph = nx.gnm_random_graph(len(genes), n_edges, seed=seed)
    random_graph = nx.relabel_nodes(random_graph, {i: g for i, g in enumerate(genes)})

    real_original_weights = {}
    for u, v in random_graph.edges():
        random_graph[u][v]['weight'] = 1.0
        real_original_weights[(u, v)] = 1.0

    return random_graph, real_original_weights
