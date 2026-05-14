import os
import pandas as pd
from scipy.stats import pearsonr, norm
from collections import defaultdict
from tqdm import tqdm
from multiprocessing import Pool
import shutil

# ====== Configuration ======

disease_name = "Paclitaxel"
dataset_name = "GSE66305"
organ = "Breast"

base_dir = f"/proj/c.zihao/work1/2drugs/{disease_name}/{dataset_name}"
output_path = f"{base_dir}/SSN"
test_filename = f"{base_dir}/data/exprSet_filtered.csv"

if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.makedirs(output_path)

# ====== Step 1: Load and Filter GTEx Data ======
Normal = pd.read_csv("/proj/c.zihao/work1/0ref/GTEx/combined_expr_df.csv", index_col=0)
organ_df = Normal[Normal["organ"] == organ]
organ_df = organ_df.drop(columns=["organ"]).T
ref_df = organ_df

# Filter out low-variance genes
initial_gene_count = ref_df.shape[0]
ref_df = ref_df.loc[ref_df.std(axis=1) > 1e-5]
filtered_gene_count = ref_df.shape[0]
print(f"Filtered constant genes: {initial_gene_count - filtered_gene_count} removed, {filtered_gene_count} remaining")

# ====== Step 2: Load and Filter Test Expression Data ======
test_df = pd.read_csv(test_filename, index_col=0)
print(test_df.head())
test_df = test_df.loc[test_df.std(axis=1) > 0]

# Align by common genes
common_genes = ref_df.index.intersection(test_df.index)

# ====== Step 3: Load and Filter STRING Network ======
links = pd.read_csv("/proj/c.zihao/work1/1survival/links.csv", index_col=0)
genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
common_genes = common_genes.intersection(genes_in_links)

# Filter links and expression matrices based on common genes
links = links[
    links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
].drop_duplicates(subset=['protein1', 'protein2'])

print(f"Filtered STRING links: {links.shape}")

# Align expression data to genes in links
used_genes = pd.unique(links[['protein1', 'protein2']].values.ravel())
ref_aligned = ref_df.loc[used_genes]
test_aligned = test_df.loc[used_genes]
print(f"Number of genes used in expression: {len(used_genes)}")

# ====== Step 4: Extract final gene pairs from STRING directly ======
gene_pairs = [
    (g1, g2)
    for g1, g2 in zip(links["protein1"], links["protein2"])
    if g1 in ref_aligned.index and g2 in ref_aligned.index
]
print(f"Total STRING gene pairs with expression data: {len(gene_pairs)}")

# ====== Step 5: Compute Reference PCC ======
n = ref_aligned.shape[1]

def compute_ref_pcc(pair):
    g1, g2 = pair
    try:
        r, _ = pearsonr(ref_aligned.loc[g1], ref_aligned.loc[g2])
    except:
        r = 0
    return (pair, r)

print("Computing reference PCCs...")

pcc_ref = dict()
for pair in tqdm(gene_pairs, desc="Reference PCC"):
    pcc_ref[pair] = compute_ref_pcc(pair)[1]


# ====== Step 6: Compute Z-score and p-value for each sample ======
def compute_sample_stat(sample_name):
    sample_df = test_aligned[[sample_name]]
    combined_df = pd.concat([ref_aligned, sample_df], axis=1)
    
    z_scores, p_values, delta_pccs = {}, {}, {}

    for edge, pcc_n in pcc_ref.items():
        g1, g2 = edge
        try:
            r, _ = pearsonr(combined_df.loc[g1], combined_df.loc[g2])
        except:
            r = 0
        delta = r - pcc_n
        sigma = (1 - pcc_n ** 2) / (n - 1)
        if sigma == 0:
            continue
        z = delta / sigma
        p = 2 * (1 - norm.cdf(abs(z)))  # Two-tailed

        z_scores[edge] = round(z, 3)
        delta_pccs[edge] = round(delta, 3)
        p_values[edge] = p

    return sample_name, z_scores, delta_pccs, p_values

print("Computing stats for test samples...")
results = []
for sample in tqdm(test_aligned.columns, desc="Samples"):
    results.append(compute_sample_stat(sample))


# ====== Step 7: Save Results ======
z_scores_dict = defaultdict(dict)
delta_pcc_dict = defaultdict(dict)
p_value_dict = defaultdict(dict)

for sample_name, z_scores, delta_pccs, p_values in results:
    for edge, z in z_scores.items():
        z_scores_dict[edge][sample_name] = z
    for edge, d in delta_pccs.items():
        delta_pcc_dict[edge][sample_name] = d
    for edge, p in p_values.items():
        p_value_dict[edge][sample_name] = p

# Convert to DataFrame
z_df = pd.DataFrame.from_dict(z_scores_dict, orient='index')
z_df.index = pd.MultiIndex.from_tuples(z_df.index, names=["Gene1", "Gene2"])

delta_df = pd.DataFrame.from_dict(delta_pcc_dict, orient='index')
delta_df.index = pd.MultiIndex.from_tuples(delta_df.index, names=["Gene1", "Gene2"])

pval_df = pd.DataFrame.from_dict(p_value_dict, orient='index')
pval_df.index = pd.MultiIndex.from_tuples(pval_df.index, names=["Gene1", "Gene2"])


z_df.to_csv(f"{output_path}/Zscore.csv")
delta_df.to_csv(f"{output_path}/delta.csv")
pval_df.to_csv(f"{output_path}/pvalue.csv")

print("ALL samples for SSN complete")
