"""
Multi-omic RWR variants for sample-specific PPI networks.

Primary model used in the 4coupled workflow:
  - MOSSN_Direct: EXP full PPI backbone + same-gene cross-layer edges to MET/CNV
    nodes, with optional dynamic concordance weighting on those cross-layer edges.

Comparative variants:
  - MOSSN_DirectNoDyn: MOSSN_Direct with fixed cross-layer edges (no concordance weighting)
  - MOSSN_Restart: restart-level coupling from non-EXP layers into EXP
  - MOSSN_NoCross: late fusion of independent single-layer runs
  - MOSSN_MultiLayer: one RWR on a full multilayer graph

Single-omic and partial-omic baselines:
  - MOSSN_EXP, MOSSN_MET, MOSSN_CNV
"""

import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix


def _norm_id(x):
    return x.replace("-", "_")[:12]


def _iqr_norm(series):
    median = series.median()
    iqr = series.quantile(0.75) - series.quantile(0.25)
    if iqr == 0:
        iqr = 1e-6
    return (series - median) / iqr


def prepare_data_MOSSN_Coupled(links, omic_data):
    """
    Build shared PPI graph on gene space = intersection of all active omics + PPI.
    Returns G, real_original_weights, omic_data (filtered & aligned).
    """
    omic_data = {
        k: v.copy().rename(columns=_norm_id)
        for k, v in omic_data.items()
    }

    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = pd.Index(genes_in_links)
    for v in omic_data.values():
        common_genes = common_genes.intersection(v.index)

    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])

    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    common_genes = common_genes.intersection(used_genes)

    common_samples = None
    for v in omic_data.values():
        s = v.columns
        common_samples = s if common_samples is None else common_samples.intersection(s)

    omic_data = {
        k: v.loc[list(common_genes), list(common_samples)].sort_index(axis=0).sort_index(axis=1)
        for k, v in omic_data.items()
    }

    G = nx.Graph()
    real_original_weights = {}
    for _, row in link_filtered.iterrows():
        u, v = row['protein1'], row['protein2']
        G.add_edge(u, v, weight=row['score'])
        real_original_weights[(u, v)] = row['score']

    return G, real_original_weights, omic_data


def _build_transition(G, real_original_weights, z_norm, alpha_mod=1.0, gamma=2.0):
    """Build omic-modulated sparse transition matrix."""
    nodes = list(G.nodes)
    n = len(nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}

    rows, cols, vals = [], [], []
    for u, v in G.edges():
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        zu = z_norm.get(u, 0)
        zv = z_norm.get(v, 0)
        sig = 1 / (1 + np.exp(-alpha_mod * (zu + zv)))
        w = max(base * (1 + gamma * (sig - 0.5)), 0)
        i, j = idx[u], idx[v]
        rows += [i, j]; cols += [j, i]; vals += [w, w]

    T = csr_matrix((vals, (rows, cols)), shape=(n, n))
    rs = np.array(T.sum(axis=1)).flatten()
    rs[rs == 0] = 1
    return T.multiply(1 / rs[:, None]), nodes, idx


def _run_rwr(T, P0, rwr_alpha=0.3, max_iter=50, tol=1e-4):
    P = P0.copy()
    for _ in range(max_iter):
        P_new = rwr_alpha * P0 + (1 - rwr_alpha) * T.T @ P
        if np.linalg.norm(P - P_new) < tol:
            return P_new
        P = P_new
    return P


def _rank_norm(P):
    ranks = np.argsort(np.argsort(P))
    return ranks / max(len(P) - 1, 1)


def _layer_rwr(sample_id, omic_data, omic_name, G, real_original_weights, rwr_alpha=0.3):
    """
    Run independent MOSSN_noPrior RWR on one omic layer.
    Returns importance vector (array, aligned to G.nodes order) and z_norm (dict).
    """
    nodes = list(G.nodes)
    data = omic_data[omic_name]
    if sample_id not in data.columns:
        return np.ones(len(nodes)) / len(nodes), {}

    z_norm = _iqr_norm(data[sample_id]).to_dict()
    T, nodes, idx = _build_transition(G, real_original_weights, z_norm)

    q90 = data[sample_id].quantile(0.9)
    seed_genes = set(data[sample_id][data[sample_id] > q90].index).intersection(G.nodes)
    P0 = np.array([1.0 if nd in seed_genes else 0.0 for nd in nodes])
    if P0.sum() == 0:
        P0[:] = 1.0
    P0 /= P0.sum()

    P = _run_rwr(T, P0, rwr_alpha)
    return _rank_norm(P), z_norm


def _edge_df(sample_id, G, real_original_weights, z_exp, importance, nodes, idx):
    """Build output DataFrame from EXP layer results."""
    rows = []
    for u, v in G.edges():
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        zu, zv = z_exp.get(u, 0), z_exp.get(v, 0)
        sig = 1 / (1 + np.exp(-(zu + zv)))
        adj_w = max(base * (1 + 2 * (sig - 0.5)), 0)
        fw = adj_w * (importance[idx[u]] + importance[idx[v]])
        rows.append({"Sample": sample_id, "Node1": u, "Node2": v,
                     "BaseWeight": base, "FinalWeight": round(fw, 5)})
    return pd.DataFrame(rows)


def _build_sparse_from_adjusted_graph(G):
    """Convert a graph with adjusted_weight on edges into a row-normalized sparse matrix."""
    nodes = list(G.nodes)
    n = len(nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}

    rows, cols, vals = [], [], []
    for u, v, data in G.edges(data=True):
        w = data['adjusted_weight']
        i, j = idx[u], idx[v]
        rows += [i, j]
        cols += [j, i]
        vals += [w, w]

    T = csr_matrix((vals, (rows, cols)), shape=(n, n))
    rs = np.array(T.sum(axis=1)).flatten()
    rs[rs == 0] = 1
    return T.multiply(1 / rs[:, None]), nodes, idx


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def MOSSN_Coupled_single_sample(sample_id, G, real_original_weights, omic_data,
                                coupled_omics, dynamic_cross=True,
                                beta=0.7, rwr_alpha=0.3):
    """
    Restart-coupled comparative variant.

    This is the implementation used by MOSSN_Restart. Non-EXP omics are first
    propagated independently, then injected into the EXP restart distribution.

    Parameters
    ----------
    coupled_omics : list of str, e.g. ['MET', 'CNV']
        Non-EXP omics whose importance feeds into EXP restart distribution.
    dynamic_cross : bool
        True  -> coupling_l[g] = importance_l[g] * sigmoid(z_l[g] * z_EXP[g])
        False -> coupling_l[g] = importance_l[g]  (fixed, no concordance scaling)
    beta : float
        Weight of EXP-native seed vs inter-layer signal in restart distribution.
        beta=1 -> pure EXP (equivalent to MOSSN_EXP).
    """
    nodes = list(G.nodes)
    n = len(nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}

    # ── Step 1: independent RWR on each non-EXP layer ────────────────────────
    layer_importance = {}
    layer_znorm = {}
    for omic in coupled_omics:
        if omic in omic_data:
            imp, znorm = _layer_rwr(sample_id, omic_data, omic, G,
                                    real_original_weights, rwr_alpha)
            layer_importance[omic] = imp
            layer_znorm[omic] = znorm

    # ── Step 2: EXP IQR normalization ────────────────────────────────────────
    exp_data = omic_data["EXP"]
    if sample_id not in exp_data.columns:
        return pd.DataFrame()
    z_exp = _iqr_norm(exp_data[sample_id]).to_dict()

    # ── Step 3: concordance-weighted coupling vectors ─────────────────────────
    P_ext = np.zeros(n)
    n_coupled = len(layer_importance)
    for omic, imp in layer_importance.items():
        if dynamic_cross:
            z_other = layer_znorm[omic]
            concordance = np.array([
                1 / (1 + np.exp(-(z_other.get(nd, 0) * z_exp.get(nd, 0))))
                for nd in nodes
            ])
            coupling = imp * concordance
        else:
            coupling = imp.copy()

        s = coupling.sum()
        if s > 0:
            coupling /= s
        P_ext += coupling

    if n_coupled > 0:
        P_ext /= n_coupled

    # ── Step 4: EXP layer RWR with augmented restart ──────────────────────────
    T_exp, _, _ = _build_transition(G, real_original_weights, z_exp)

    q90 = exp_data[sample_id].quantile(0.9)
    seed_genes = set(exp_data[sample_id][exp_data[sample_id] > q90].index).intersection(G.nodes)
    P0_exp = np.array([1.0 if nd in seed_genes else 0.0 for nd in nodes])
    if P0_exp.sum() == 0:
        P0_exp[:] = 1.0
    P0_exp /= P0_exp.sum()

    if n_coupled > 0 and P_ext.sum() > 0:
        P_ext /= P_ext.sum()
        P0_combined = beta * P0_exp + (1 - beta) * P_ext
    else:
        P0_combined = P0_exp

    P0_combined /= P0_combined.sum()
    P_final = _run_rwr(T_exp, P0_combined, rwr_alpha)
    importance = _rank_norm(P_final)

    return _edge_df(sample_id, G, real_original_weights, z_exp, importance, nodes, idx)


def prepare_data_MOSSN_Direct(links, omic_data, direct_omics):
    """
    Build a mixed graph:
      - EXP layer: full PPI network  (nodes = gene names)
      - MET/CNV layers: isolated nodes only  (nodes = 'MET:gene', 'CNV:gene')
      - Cross-layer edges: gene <-> 'OMIC:gene'  (base weight = 1.0)

    Returns G_direct, weights_direct, omic_data (filtered & aligned), exp_genes (list).
    """
    # filter & align omics (gene space = ∩ all omics ∩ PPI)
    all_omics = {"EXP": omic_data["EXP"]}
    for o in direct_omics:
        if o in omic_data:
            all_omics[o] = omic_data[o]

    all_omics = {
        k: v.copy().rename(columns=_norm_id)
        for k, v in all_omics.items()
    }

    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = pd.Index(genes_in_links)
    for v in all_omics.values():
        common_genes = common_genes.intersection(v.index)

    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])

    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    common_genes = common_genes.intersection(used_genes)

    common_samples = None
    for v in all_omics.values():
        s = v.columns
        common_samples = s if common_samples is None else common_samples.intersection(s)

    all_omics = {
        k: v.loc[list(common_genes), list(common_samples)].sort_index(axis=0).sort_index(axis=1)
        for k, v in all_omics.items()
    }

    G = nx.Graph()
    weights = {}

    # EXP intra-layer PPI edges
    for _, row in link_filtered.iterrows():
        u, v = row['protein1'], row['protein2']
        G.add_edge(u, v, weight=row['score'])
        weights[(u, v)] = row['score']

    # isolated MET/CNV nodes + cross-layer edges
    for omic in direct_omics:
        if omic not in all_omics:
            continue
        for gene in common_genes:
            omic_node = f"{omic}:{gene}"
            G.add_node(omic_node)
            G.add_edge(gene, omic_node, weight=1.0)
            weights[(gene, omic_node)] = 1.0

    return G, weights, all_omics, list(common_genes)


def prepare_data_MOSSN_Multilayer(links, omic_data, active_omics):
    """
    Build a full multilayer graph:
      - Each omic has its own PPI layer: 'OMIC:gene' nodes with intra-layer PPI edges
      - Same gene is connected only between EXP and each non-EXP omic layer:
        'EXP:gene' <-> 'MET:gene', 'EXP:gene' <-> 'CNV:gene'

    Returns G_multi, weights_multi, omic_data (filtered & aligned), exp_genes (list).
    """
    all_omics = {
        k: omic_data[k]
        for k in active_omics
        if k in omic_data
    }
    if "EXP" not in all_omics:
        raise ValueError("EXP must be included in active_omics for multilayer RWR.")

    all_omics = {
        k: v.copy().rename(columns=_norm_id)
        for k, v in all_omics.items()
    }

    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = pd.Index(genes_in_links)
    for v in all_omics.values():
        common_genes = common_genes.intersection(v.index)

    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
    ].drop_duplicates(subset=['protein1', 'protein2'])

    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    common_genes = common_genes.intersection(used_genes)

    common_samples = None
    for v in all_omics.values():
        s = v.columns
        common_samples = s if common_samples is None else common_samples.intersection(s)

    all_omics = {
        k: v.loc[list(common_genes), list(common_samples)].sort_index(axis=0).sort_index(axis=1)
        for k, v in all_omics.items()
    }

    G = nx.Graph()
    weights = {}
    active_omics = list(all_omics.keys())

    # Intra-layer PPI edges for every active omic.
    for omic in active_omics:
        for _, row in link_filtered.iterrows():
            u = f"{omic}:{row['protein1']}"
            v = f"{omic}:{row['protein2']}"
            G.add_edge(u, v, weight=row['score'], edge_type="intra", omic=omic)
            weights[(u, v)] = row['score']

    # Inter-layer same-gene edges only between EXP and each non-EXP layer.
    for gene in common_genes:
        exp_node = f"EXP:{gene}"
        for omic in active_omics:
            if omic == "EXP":
                continue
            omic_node = f"{omic}:{gene}"
            G.add_edge(exp_node, omic_node, weight=1.0, edge_type="inter", gene=gene)
            weights[(exp_node, omic_node)] = 1.0

    return G, weights, all_omics, list(common_genes)


def MOSSN_Direct_single_sample(sample_id, G, real_original_weights, omic_data,
                               exp_genes, direct_omics,
                               dynamic_cross=True, rwr_alpha=0.3):
    """
    Direct cross-layer RWR.

    This is the primary model used in the 4coupled workflow.
      - EXP intra-layer edges modulated by EXP z-scores (sigmoid)
      - Cross-layer edges gene <-> OMIC:gene modulated by concordance
          dynamic: w = sigmoid(z_omic[g] * z_EXP[g])
          fixed  : w = 1.0  (base weight unchanged)
      - MET/CNV nodes are isolated (no intra-layer edges)
      - Seed nodes from all active omics
      - FinalWeight computed from EXP-layer node importance only
    """
    alpha_mod = 1.0
    gamma     = 2.0

    # ── IQR normalization per omic ────────────────────────────────────────────
    z_norms = {}
    for omic, data in omic_data.items():
        if sample_id in data.columns:
            z_norms[omic] = _iqr_norm(data[sample_id]).to_dict()

    z_exp = z_norms.get("EXP", {})

    # ── edge weight modulation ────────────────────────────────────────────────
    for u, v, data in G.edges(data=True):
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))

        u_is_exp = ":" not in u
        v_is_exp = ":" not in v

        if u_is_exp and v_is_exp:
            # EXP intra-layer: sigmoid of sum of EXP z-scores
            zu = z_exp.get(u, 0)
            zv = z_exp.get(v, 0)
            sig = 1 / (1 + np.exp(-alpha_mod * (zu + zv)))
            data['adjusted_weight'] = max(base * (1 + gamma * (sig - 0.5)), 0)

        else:
            # cross-layer edge: gene (EXP) <-> OMIC:gene
            gene     = u if u_is_exp else v
            omic_nd  = v if u_is_exp else u
            omic     = omic_nd.split(":")[0]
            z_other  = z_norms.get(omic, {})

            if dynamic_cross:
                zo  = z_other.get(gene, 0)
                ze  = z_exp.get(gene, 0)
                sig = 1 / (1 + np.exp(-(zo * ze)))
                data['adjusted_weight'] = max(base * (1 + gamma * (sig - 0.5)), 0)
            else:
                data['adjusted_weight'] = base

    # ── sparse transition matrix ──────────────────────────────────────────────
    nodes = list(G.nodes)
    n     = len(nodes)
    idx   = {nd: i for i, nd in enumerate(nodes)}

    rows, cols, vals = [], [], []
    for u, v, data in G.edges(data=True):
        w    = data['adjusted_weight']
        i, j = idx[u], idx[v]
        rows += [i, j]; cols += [j, i]; vals += [w, w]

    T  = csr_matrix((vals, (rows, cols)), shape=(n, n))
    rs = np.array(T.sum(axis=1)).flatten()
    rs[rs == 0] = 1
    T  = T.multiply(1 / rs[:, None])

    # ── seed nodes from all active omics ─────────────────────────────────────
    seed_nodes = set()
    # EXP: top-90% expression genes
    exp_data = omic_data["EXP"]
    if sample_id in exp_data.columns:
        q90 = exp_data[sample_id].quantile(0.9)
        seed_nodes.update(g for g in exp_data[sample_id][exp_data[sample_id] > q90].index
                          if g in G.nodes)
    # MET/CNV: top-90% of their values -> OMIC:gene nodes
    for omic in direct_omics:
        if omic not in omic_data or sample_id not in omic_data[omic].columns:
            continue
        s   = omic_data[omic][sample_id]
        q90 = s.quantile(0.9)
        seed_nodes.update(f"{omic}:{g}" for g in s[s > q90].index
                          if f"{omic}:{g}" in G.nodes)

    P0 = np.array([1.0 if nd in seed_nodes else 0.0 for nd in nodes])
    if P0.sum() == 0:
        P0[:] = 1.0
    P0 /= P0.sum()

    # ── RWR ──────────────────────────────────────────────────────────────────
    P = _run_rwr(T, P0, rwr_alpha)

    # rank-normalize using only EXP-layer nodes
    exp_idx   = [idx[g] for g in exp_genes if g in idx]
    P_exp_sub = P[exp_idx]
    ranks_sub = np.argsort(np.argsort(P_exp_sub))
    norm_sub  = ranks_sub / max(len(ranks_sub) - 1, 1)
    # map back to full node importance
    importance = np.zeros(n)
    for rank_i, node_i in enumerate(exp_idx):
        importance[node_i] = norm_sub[rank_i]

    # ── output: EXP intra-layer edges only ───────────────────────────────────
    rows_out = []
    for u, v, data in G.edges(data=True):
        if ":" in u or ":" in v:
            continue   # skip cross-layer edges in output
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        fw   = data['adjusted_weight'] * (importance[idx[u]] + importance[idx[v]])
        rows_out.append({"Sample": sample_id, "Node1": u, "Node2": v,
                         "BaseWeight": base, "FinalWeight": round(fw, 5)})
    return pd.DataFrame(rows_out)


def MOSSN_Multilayer_single_sample(sample_id, G, real_original_weights, omic_data,
                                   exp_genes, active_omics,
                                   dynamic_cross=True, rwr_alpha=0.3):
    """
    Full multilayer RWR:
      - Every active omic has its own PPI layer
      - Same genes are linked across omic layers
      - One RWR runs on the entire multilayer graph
      - Output keeps only EXP-layer intra-layer edges
    """
    alpha_mod = 1.0
    gamma = 2.0

    z_norms = {}
    for omic, data in omic_data.items():
        if sample_id in data.columns:
            z_norms[omic] = _iqr_norm(data[sample_id]).to_dict()

    for u, v, data in G.edges(data=True):
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        u_omic, u_gene = u.split(":", 1)
        v_omic, v_gene = v.split(":", 1)

        if data.get("edge_type") == "intra":
            z_u = z_norms.get(u_omic, {}).get(u_gene, 0)
            z_v = z_norms.get(v_omic, {}).get(v_gene, 0)
            sig = 1 / (1 + np.exp(-alpha_mod * (z_u + z_v)))
            data['adjusted_weight'] = max(base * (1 + gamma * (sig - 0.5)), 0)
        else:
            z_u = z_norms.get(u_omic, {}).get(u_gene, 0)
            z_v = z_norms.get(v_omic, {}).get(v_gene, 0)
            if dynamic_cross:
                sig = 1 / (1 + np.exp(-(z_u * z_v)))
                data['adjusted_weight'] = max(base * (1 + gamma * (sig - 0.5)), 0)
            else:
                data['adjusted_weight'] = base

    T, nodes, idx = _build_sparse_from_adjusted_graph(G)

    seed_nodes = set()
    for omic in active_omics:
        if omic not in omic_data or sample_id not in omic_data[omic].columns:
            continue
        s = omic_data[omic][sample_id]
        q90 = s.quantile(0.9)
        seed_nodes.update(
            f"{omic}:{gene}"
            for gene in s[s > q90].index
            if f"{omic}:{gene}" in G.nodes
        )

    P0 = np.array([1.0 if nd in seed_nodes else 0.0 for nd in nodes])
    if P0.sum() == 0:
        P0[:] = 1.0
    P0 /= P0.sum()

    P = _run_rwr(T, P0, rwr_alpha)

    exp_nodes = [f"EXP:{gene}" for gene in exp_genes if f"EXP:{gene}" in idx]
    exp_idx = [idx[node] for node in exp_nodes]
    P_exp_sub = P[exp_idx]
    ranks_sub = np.argsort(np.argsort(P_exp_sub))
    norm_sub = ranks_sub / max(len(ranks_sub) - 1, 1)

    importance = np.zeros(len(nodes))
    for rank_i, node_i in enumerate(exp_idx):
        importance[node_i] = norm_sub[rank_i]

    z_exp = z_norms.get("EXP", {})
    rows_out = []
    for u, v, data in G.edges(data=True):
        if data.get("edge_type") != "intra":
            continue
        if not (u.startswith("EXP:") and v.startswith("EXP:")):
            continue

        gene_u = u.split(":", 1)[1]
        gene_v = v.split(":", 1)[1]
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        zu = z_exp.get(gene_u, 0)
        zv = z_exp.get(gene_v, 0)
        sig = 1 / (1 + np.exp(-(zu + zv)))
        adj_w = max(base * (1 + gamma * (sig - 0.5)), 0)
        fw = adj_w * (importance[idx[u]] + importance[idx[v]])
        rows_out.append({
            "Sample": sample_id,
            "Node1": gene_u,
            "Node2": gene_v,
            "BaseWeight": base,
            "FinalWeight": round(fw, 5),
        })
    return pd.DataFrame(rows_out)


def MOSSN_LateFusion_single_sample(sample_id, G, real_original_weights, omic_data,
                                   active_omics, rwr_alpha=0.3):
    """
    Late fusion: run independent RWR on each omic layer, average FinalWeights.
    Used for the MOSSN_NoCross comparative variant.
    """
    nodes = list(G.nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}

    all_fw = []
    for omic in active_omics:
        if omic not in omic_data:
            continue
        imp, z_norm = _layer_rwr(sample_id, omic_data, omic, G,
                                  real_original_weights, rwr_alpha)
        rows = []
        for u, v in G.edges():
            base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
            zu, zv = z_norm.get(u, 0), z_norm.get(v, 0)
            sig = 1 / (1 + np.exp(-(zu + zv)))
            adj_w = max(base * (1 + 2 * (sig - 0.5)), 0)
            fw = adj_w * (imp[idx[u]] + imp[idx[v]])
            rows.append(fw)
        all_fw.append(rows)

    if not all_fw:
        return pd.DataFrame()

    avg_fw = np.mean(all_fw, axis=0)
    out_rows = []
    for i, (u, v) in enumerate(G.edges()):
        base = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        out_rows.append({"Sample": sample_id, "Node1": u, "Node2": v,
                         "BaseWeight": base, "FinalWeight": round(avg_fw[i], 5)})
    return pd.DataFrame(out_rows)
