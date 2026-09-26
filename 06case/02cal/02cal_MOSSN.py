import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append("/proj/c.zihao/work1/function/ablation")
import MOSSN_func as MOSSN

base_dir = "/proj/c.zihao/work1/06case"
save_path = os.path.join(base_dir, "04results", "MOSSN")
expr_file = os.path.join(base_dir, "01data", "exprSet_filtered.csv")
link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)
os.chdir(save_path)

links = pd.read_csv(link_file, index_col=0)
expression_data = pd.read_csv(expr_file, index_col=0)

G, uniform_original_weights, expression_data = MOSSN.MOSSN_prepare_data(
    links=links,
    expression_data=expression_data,
)
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
print(f"Samples: {expression_data.shape[1]}")

lam = 2.0
rwr_alpha = 0.3
seed_quantile = 0.9

for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
    edge_weights_df = MOSSN.MOSSN_single_sample(
        sample_id=sample_id,
        G=G,
        uniform_original_weights=uniform_original_weights,
        expression_data=expression_data,
        lam=lam,
        rwr_alpha=rwr_alpha,
        seed_quantile=seed_quantile,
    )
    edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
