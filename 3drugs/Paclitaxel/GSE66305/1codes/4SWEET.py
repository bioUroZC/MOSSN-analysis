import os
import time
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import shutil

# ------------------------- Set input paths -------------------------

# === Set dataset name here ===
disease_name = "Paclitaxel"
dataset_name = "GSE66305"

# === Setup Working Directory and Load Global Data ===
base_dir = f"/proj/c.zihao/work1/2drugs/{disease_name}/{dataset_name}"
file_e = f"{base_dir}/data/exprSet_filtered.csv"
save_path = f"{base_dir}/SWEET"
link_file = '/proj/c.zihao/work1/1survival/links.csv'

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)

# ------------------------- Step 1: Load expression matrix -------------------------
print("[INFO] Loading expression matrix...")
expr_df = pd.read_csv(file_e, index_col=0)
expr_df.index = expr_df.index.astype(str)
expr_df.columns = expr_df.columns.astype(str)

samples = expr_df.columns.tolist()
genes = expr_df.index.tolist()
print(f"[INFO] Number of samples: {len(samples)}, number of genes: {len(genes)}")

# ------------------------- Step 2: Load gene pairs of interest from link.csv -------------------------
gene_pairs_set = set()
if os.path.exists(link_file):
    link_df = pd.read_csv(link_file, index_col=0)
    for _, row in link_df.iterrows():
        g1 = str(row['protein1']).strip()
        g2 = str(row['protein2']).strip()
        gene_pairs_set.add((g1, g2))
    print(f"[INFO] Loaded {len(gene_pairs_set)} gene pairs from link.csv")
else:
    gene_pairs_set = None
    print("[INFO] link.csv not found, all gene-gene pairs will be included")

# ------------------------- Step 3: Compute sample-sample PCC matrix -------------------------
def compute_pcc_matrix(df):
    n = df.shape[1]
    pcc_matrix = pd.DataFrame(np.zeros((n, n)), index=df.columns, columns=df.columns)
    for i in range(n):
        print(f"[INFO] Computing PCC for sample: {df.columns[i]} ({i+1}/{n})")
        for j in range(i, n):
            r, _ = pearsonr(df.iloc[:, i], df.iloc[:, j])
            pcc_matrix.iloc[i, j] = r
            pcc_matrix.iloc[j, i] = r
    return pcc_matrix

print("[INFO] Computing sample PCC matrix...")
t0 = time.time()
pcc_matrix = compute_pcc_matrix(expr_df)
print(f"[INFO] PCC matrix computed in {time.time() - t0:.2f} seconds.")

# ------------------------- Step 4: Compute sample weights -------------------------
def compute_sample_weights(pcc_matrix, k=0.1, x=0.01):
    n = pcc_matrix.shape[0]
    mean_pcc = (pcc_matrix.sum(axis=1) - 1) / (n - 1)
    rmin = mean_pcc.min()
    rmax = mean_pcc.max()
    weights = (mean_pcc - rmin + x) / (rmax - rmin + x)
    weights = weights * k * n
    return weights

weights = compute_sample_weights(pcc_matrix)
weights_dict = weights.to_dict()

weights_df = weights.reset_index()
weights_df.columns = ['patient', 'sample_weight']
weights_df.to_csv(os.path.join(save_path, 'weights.txt'), sep="\t", index=False)

# ------------------------- Step 5: Remove genes with all-zero expression -------------------------
value = expr_df.values.astype(float)
gene_arr = np.array(genes)

zero_gene_idx = np.where(np.sum(value, axis=1) == 0)[0]
if len(zero_gene_idx) > 0:
    print('[INFO] Removing genes with all-zero expression: ' + ','.join(gene_arr[zero_gene_idx]))
    value = np.delete(value, zero_gene_idx, axis=0)
    gene_arr = np.delete(gene_arr, zero_gene_idx)

# ------------------------- Step 6: Compute global gene-gene PCC -------------------------
agg = np.corrcoef(value)

# ------------------------- Step 7: Build sample-specific networks (filtered by link.csv) -------------------------
all_scores = []
for idx, sample in enumerate(samples):
    if sample not in weights_dict:
        continue

    expr_aug = np.c_[value, value[:, idx]]
    pcc_s = np.corrcoef(expr_aug)
    w = weights_dict[sample]
    weighted_pcc = w * (pcc_s - agg) + agg

    edges = []
    for i in range(len(gene_arr)):
        for j in range(i + 1, len(gene_arr)):
            g1, g2 = gene_arr[i], gene_arr[j]
            if gene_pairs_set is not None and (g1, g2) not in gene_pairs_set and (g2, g1) not in gene_pairs_set:
                continue
            score = weighted_pcc[i, j]
            edges.append([g1, g2, score])
            all_scores.append(score)

    df_out = pd.DataFrame(edges, columns=["gene1", "gene2", "raw_edge_score"])
    df_out.to_csv(f"{save_path}/{sample}.txt", sep="\t", index=False)

print("[INFO] All sample-specific networks saved.")

# ------------------------- Step 8: Calculate and save z-score for all edges -------------------------
all_scores = np.array(all_scores)
vmean = np.mean(all_scores)
vstd = np.std(all_scores)

with open(os.path.join(save_path, "mean_std.txt"), mode='w') as wline:
    wline.write(f"mean\t{vmean}\nstd\t{vstd}\n")

for sample in samples:
    net_file = os.path.join(save_path, f"{sample}.txt")
    if not os.path.exists(net_file):
        continue

    df_net = pd.read_csv(net_file, sep="\t")
    df_net['z_score'] = (df_net['raw_edge_score'] - vmean) / vstd
    df_net[['gene1', 'gene2', 'z_score']].to_csv(
        os.path.join(save_path, f"{sample}_zscore.txt"),
        sep="\t", index=False
    )

print("[INFO] All sameple for SWEET saved.")
