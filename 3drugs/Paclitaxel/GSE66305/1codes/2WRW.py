import os
import time
import glob
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from scipy.sparse import csr_matrix
import shutil

# === Set dataset name here ===

disease_name = "Paclitaxel"
dataset_name = "GSE66305"

# === Setup Working Directory and Load Global Data ===
base_dir = f"/proj/c.zihao/work1/2drugs/{disease_name}/{dataset_name}"
save_path = f"{base_dir}/WRW"

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)

os.chdir(save_path)

# Load seed nodes (SS markers)
seed_nodes_df = pd.read_csv("/proj/c.zihao/work1/1survival/SSmarkers.csv", index_col=0)
print(seed_nodes_df.head())
print("Seed Data Loaded.")

# Load links (protein interactions)
links = pd.read_csv("/proj/c.zihao/work1/1survival/links.csv", index_col=0)
print(links.head())
print(links.shape)
print("Links Data Loaded.")

# Load expression data
expression_data = pd.read_csv(f"{base_dir}/data/exprSet_filtered.csv", index_col=0)
print(expression_data.iloc[1:5, 1:5])
print("Successfully read exprSet_filtered.csv")

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
start_time = time.time()
G = nx.Graph()
for _, row in links.iterrows():
    G.add_edge(row['protein1'], row['protein2'], weight=row['score'])

# Store real original weights for each edge
real_original_weights = {
    (row['protein1'], row['protein2']): row['score']
    for _, row in links.iterrows()
}
step_time = time.time() - start_time
print("Network created with nodes and edges.")
print(f"Number of nodes: {G.number_of_nodes()}, Number of edges: {G.number_of_edges()}")
print(f"Step 1 (Construct Protein Interaction Network ) Time: {step_time:.4f} seconds")

# === Process Each Sample ===
# Iterate through each sample with a progress bar
for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
    print(f"\nProcessing Sample: {sample_id}")

    # Reset edge weights to real original values
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):  # Check to avoid missing edges
            G[u][v]['weight'] = weight

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
        modulation_factor = 1 + gamma * (sigmoid - 0.5)
        modulation_factor = max(modulation_factor, 0)
        data['weight'] = real_original_weights.get((u, v), real_original_weights.get((v, u), 0)) * modulation_factor
        data['modulation_factor'] = modulation_factor
        

    step_time = time.time() - start_time
    print(f"Step 3 (Modulation-Adjusted Edge Weights) Time: {step_time:.4f} seconds")

    # === Step 4: Build Transition Matrix ===
    start_time = time.time()
    #seed_nodes = {gene for gene, expr in sample_expression.items() if expr > sample_mean}

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
        weight = G[u][v]['weight']
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
    rwr_alpha = 0.3   # Reset probability
    P0 = np.array([1 if node in seed_nodes else 0 for node in nodes])
    P0 = P0 / P0.sum()  # Normalize
    P = P0.copy()
    max_iter = 50  # Limit the number of iterations
    tol = 1e-4  # Relaxed convergence tolerance

    for step in range(max_iter):
        P_new = rwr_alpha * P0 + (1 - rwr_alpha) * T_sparse.T @ P  # Sparse matrix-vector multiplication
        if np.linalg.norm(P - P_new) < tol:  # Dense norm for convergence check
            print(f"Converged at Step {step + 1}")
            break
        P = P_new

    # Sort P and get the ranks
    ranks = np.argsort(np.argsort(P))
    # Normalize using ranks to ensure an even distribution between 0 and 1
    P_normalized = ranks / (len(P) - 1)

    step_time = time.time() - start_time
    print(f"Step 5 (Random Walk Calculation) Time: {step_time:.4f} seconds")

    # Save Node Scores for Current Sample
    start_time = time.time()
    node_scores = pd.DataFrame({
        "Node": nodes,
        "Score": P
    }).sort_values(by="Score", ascending=False)
    node_scores["Sample"] = sample_id
    node_scores.to_csv(f"{sample_id}_nodes.csv", index=False)
    step_time = time.time() - start_time
    print(f"Step 6 (Save Node Scores) Time: {step_time:.4f} seconds")

    # Save updated edge weights after random walk
    start_time = time.time()
    edges_data = []

    for u, v, data in G.edges(data=True):
        original_weight = data['weight']  # Dynamically adjusted edge weight
        importance_u = P_normalized[node_to_index[u]]  # Probability of node u from random walk
        importance_v = P_normalized[node_to_index[v]]  # Probability of node v from random walk

        # Get the expression values of Node1 and Node2
        expression_u = expression_data.loc[u, sample_id] if u in expression_data.index else np.nan
        expression_v = expression_data.loc[v, sample_id] if v in expression_data.index else np.nan

        # Calculate the new weight using a chosen formula (e.g., sum of probabilities)
        updated_weight = original_weight * (importance_u + importance_v)

        # Retrieve the score (real original weight)
        score = real_original_weights.get((u, v), real_original_weights.get((v, u), 0))  # Handle undirected edges
        
        modulation_factor = G[u][v].get('modulation_factor', np.nan)

        # Store the edge information
        edges_data.append({
            "Sample": sample_id,
            "Node1": u,
            "Node2": v,
            "Expression_Node1": expression_u, 
            "Expression_Node2": expression_v, 
            "Importance_Node1": importance_u, 
            "Importance_Node2": importance_v,  
            "PPIWeight": score,
            "ModulationFactor": modulation_factor,  
            "UpdatedWeight": original_weight,
            "FinalWeight": updated_weight     
        })

    edge_weights_df = pd.DataFrame(edges_data)
    print(edge_weights_df.head())
    print(edge_weights_df.shape)
    edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)

    step_time = time.time() - start_time
    print(f"Step 7 (Save Edge Weights) Time: {step_time:.4f} seconds")

print("All samples for WRW completed")
