import time
import glob
import shutil
import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix


def prepare_data_RWR(seed_nodes_df, links, expression_data):
    expression_data = expression_data.loc[expression_data.std(axis=1) > 0]
    # Now get common genes with updated expression_data
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    # Filter links based on updated gene list
    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
        ].drop_duplicates(subset=['protein1', 'protein2'])
    links = link_filtered
    print(links.shape)
    # Filter expression matrix again (just to be safe)
    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]
    # === Construct Protein Interaction Network ===

    G = nx.Graph()
    for _, row in links.iterrows():
        G.add_edge(row['protein1'], row['protein2'], weight=row['score'])

    # Store real original weights for each edge
    real_original_weights = {
        (row['protein1'], row['protein2']): row['score']
        for _, row in links.iterrows()
    }

    return G, real_original_weights, expression_data, seed_nodes_df



def RWR_single_sample(sample_id, G, real_original_weights, expression_data, seed_nodes_df):
    print(f"\nProcessing Sample: {sample_id}")

    # Reset edge weights to real original values
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):  # Check to avoid missing edges
            G[u][v]['weight'] = weight
            G[u][v]['adjusted_weight'] = weight
            G[u][v]['correction_score'] = 1.0

    # === Step 2: Select Sample Data ===
    start_time = time.time()
    sample_expression = expression_data[sample_id]
    step_time = time.time() - start_time
    print(f"Step 2 (Select Sample Data) Time: {step_time:.4f} seconds")

    # === Step 3: Adjust Edge Weights Dynamically (Using Modulation Score and Sigmoid) ===
    start_time = time.time()

    # Step 3.1: Compute sample-specific median and IQR for expression
    sample_median = sample_expression.median()
    print(sample_median)
    sample_iqr = sample_expression.quantile(0.75) - sample_expression.quantile(0.25)
    if sample_iqr == 0:
        sample_iqr = 1e-6  # prevent division by zero

    # Step 3.2: Compute normalized expression deviation
    expr_norm = ((sample_expression - sample_median) / sample_iqr).copy()

    # Step 3.3: Adjust edge weights using modulation score and sigmoid
    alpha_mod = 1.0
    gamma = 2.0


    for u, v, data in G.edges(data=True):
        zu = expr_norm.get(u, 0)
        zv = expr_norm.get(v, 0)
        M_uv = zu + zv
        sigmoid = 1 / (1 + np.exp(-alpha_mod * M_uv))
        correction_score = 1 + gamma * (sigmoid - 0.5)
        correction_score = max(correction_score, 0)
        data['adjusted_weight'] = real_original_weights.get((u, v),
                                                            real_original_weights.get((v, u), 0)) * correction_score
        data['correction_score'] = correction_score

    step_time = time.time() - start_time
    print(f"Step 3 (Modulation-Adjusted Edge Weights) Time: {step_time:.4f} seconds")

    # === Step 4: Build Transition Matrix ===
    start_time = time.time()
    # seed_nodes = {gene for gene, expr in sample_expression.items() if expr > sample_mean}

    # Extract seed nodes from the 'symbol' column
    ssmarker_genes = set(seed_nodes_df['symbol'])
    q90 = sample_expression.quantile(0.9)
    high_expr_genes = set(sample_expression[sample_expression > q90].index)
    print(f"high_expr_gene number: {len(high_expr_genes)}")
    print("high_expr_genes:", list(high_expr_genes)[:10])
    seed_nodes = ssmarker_genes.union(high_expr_genes)

    # Ensure the seed nodes exist in the graph
    seed_nodes = seed_nodes.intersection(G.nodes)
    print(f"Number of valid seed nodes: {len(seed_nodes)}")

    nodes = list(G.nodes)
    n = len(nodes)

    # Create a mapping from node to index for faster lookup
    node_to_index = {node: i for i, node in enumerate(nodes)}

    # Prepare sparse matrix data
    row, col, edge_weights = [], [], []
    for u, v in G.edges():
        weight = G[u][v]['adjusted_weight']
        row.append(node_to_index[u])
        col.append(node_to_index[v])
        edge_weights.append(weight)

        row.append(node_to_index[v])
        col.append(node_to_index[u])
        edge_weights.append(weight)

    # Create and normalize the sparse transition matrix
    T_sparse = csr_matrix((edge_weights, (row, col)), shape=(n, n))
    row_sums = np.array(T_sparse.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    T_sparse = T_sparse.multiply(1 / row_sums[:, None])  # Normalize rows
    step_time = time.time() - start_time
    print(f"Step 4 (Build Transition Matrix) Time: {step_time:.4f} seconds")

    # === Step 5: Random Walk Calculation ===
    start_time = time.time()
    rwr_alpha = 0.3  # Reset probability
    P0 = np.array([1 if node in seed_nodes else 0 for node in nodes])

    if P0.sum() == 0:
        P0[:] = 1

    P0 = P0 / P0.sum()  # Normalize
    P = P0.copy()
    max_iter = 50  # Limit the number of iterations
    tol = 1e-4  # Relaxed convergence tolerance

    for step in range(max_iter):
        # Sparse matrix-vector multiplication
        P_new = rwr_alpha * P0 + (1 - rwr_alpha) * T_sparse.T @ P
        # Dense norm for convergence check
        if np.linalg.norm(P - P_new) < tol:
            print(f"Converged at Step {step + 1}")
            #P = P_new
            break
        P = P_new

    # Sort P and get the ranks
    ranks = np.argsort(np.argsort(P))
    # Normalize using ranks to ensure an even distribution between 0 and 1
    P_normalized = ranks / (len(P) - 1)
    pi_prob = P / P.sum()

    step_time = time.time() - start_time
    print(f"Step 5 (Random Walk Calculation) Time: {step_time:.4f} seconds")

    # Step 6:  Save Node Scores for Current Sample
    start_time = time.time()
    node_scores = pd.DataFrame({
        "Node": nodes,
        "Score": P
    }).sort_values(by="Score", ascending=False)
    node_scores["Sample"] = sample_id
    # node_scores.to_csv(f"{sample_id}_nodes.csv", index=False)
    step_time = time.time() - start_time
    print(f"Step 6 (Save Node Scores) Time: {step_time:.4f} seconds")

    # Save updated edge weights after random walk
    start_time = time.time()
    edges_data = []

    T_coo = T_sparse.tocoo()
    T_lookup = {(ri, cj): val for ri, cj, val in zip(T_coo.row, T_coo.col, T_coo.data)}

    def Tij(i, j):
        return T_lookup.get((i, j), 0.0)

    for u, v, data in G.edges(data=True):
        # Get the expression values of Node1 and Node2
        expression_u = expression_data.loc[u, sample_id] if u in expression_data.index else np.nan
        expression_v = expression_data.loc[v, sample_id] if v in expression_data.index else np.nan
        # Probability of node u from random walk
        importance_u = P_normalized[node_to_index[u]]
        # Probability of node v from random walk
        importance_v = P_normalized[node_to_index[v]]

        # Retrieve the score (real original weight)
        BaseWeight = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))
        correction_score = G[u][v].get('correction_score', np.nan)
        adjusted_weight = data['adjusted_weight']

        i = node_to_index[u]
        j = node_to_index[v]
        P_ij = T_lookup.get((i, j), 0.0)
        P_ji = T_lookup.get((j, i), 0.0)
        edge_flow = pi_prob[i] * P_ij + pi_prob[j] * P_ji  # 单纯流量

        # Store the edge information
        edges_data.append({
            "Sample": sample_id,
            "Node1": u,
            "Node2": v,
            "ExpNode1": expression_u,
            "ExpNode2": expression_v,
            "ImpNode1": importance_u,
            "ImpNode2": importance_v,
            "BaseWeight": BaseWeight,
            "CorrectionScore": correction_score,
            "AdjustedWeight": adjusted_weight,
            "Flow": edge_flow
        })

    edge_weights_df = pd.DataFrame(edges_data)

    edge_weights_df["ExpNode1"] = edge_weights_df["ExpNode1"].round(5)
    edge_weights_df["ExpNode2"] = edge_weights_df["ExpNode2"].round(5)
    edge_weights_df["ImpNode1"] = edge_weights_df["ImpNode1"].round(5)
    edge_weights_df["ImpNode2"] = edge_weights_df["ImpNode2"].round(5)
    edge_weights_df["CorrectionScore"] = edge_weights_df["CorrectionScore"].round(5)
    edge_weights_df["AdjustedWeight"] = edge_weights_df["AdjustedWeight"].round(5)

    edge_weights_df["IMPscore"] = edge_weights_df["ImpNode1"] + edge_weights_df["ImpNode2"]

    edge_weights_df["RefinedWeight"] = edge_weights_df["AdjustedWeight"] * (
            edge_weights_df["ImpNode1"] + edge_weights_df["ImpNode2"])
    edge_weights_df["RefinedWeight"] = edge_weights_df["RefinedWeight"].round(5)

    edge_weights_df["FlowWeight"] = edge_weights_df["Flow"] * edge_weights_df["AdjustedWeight"]

    edge_weights_df["Flow"] = edge_weights_df["Flow"].round(10)
    edge_weights_df["FlowWeight"] = edge_weights_df["FlowWeight"].round(10)

    edge_weights_df["IMPscoreRank"] = edge_weights_df["IMPscore"].rank(pct=True, method="average").round(5)
    edge_weights_df["RefinedWeightRank"] = edge_weights_df["RefinedWeight"].rank(pct=True, method="average").round(5)
    edge_weights_df["FlowRank"] = edge_weights_df["Flow"].rank(pct=True, method="average").round(5)
    edge_weights_df["FlowWeightRank"] = edge_weights_df["FlowWeight"].rank(pct=True, method="average").round(5)

    return edge_weights_df, node_scores