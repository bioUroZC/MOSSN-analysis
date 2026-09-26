import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import MOSSN_func as MOSSN

dataset_name = "LUAD"
base_dir = "/proj/c.zihao/work1/02analysis/01robust"
LEVELS = [10, 30, 50, 70]

lam = 2.0
rwr_alpha = 0.3
seed_quantile = 0.9

links = pd.read_csv("/proj/c.zihao/work1/01data/06string/links.csv", index_col=0)

for level in LEVELS:
    print(f"#========== {dataset_name} {level}% ==========")
    level_dir = f"{base_dir}/{level}"
    save_path = f"{level_dir}/MOSSN/{dataset_name}"
    expr_file = f"{level_dir}/data/LUAD_exprSet_half.csv"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    expression_data = pd.read_csv(expr_file, index_col=0)

    G, uniform_original_weights, expression_data = MOSSN.MOSSN_prepare_data(
        links=links,
        expression_data=expression_data
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    for sample_id in tqdm(expression_data.columns, desc=f"Processing {level}% samples", unit="sample"):
        edge_weights_df = MOSSN.MOSSN_single_sample(
            sample_id=sample_id,
            G=G,
            uniform_original_weights=uniform_original_weights,
            expression_data=expression_data,
            lam=lam,
            rwr_alpha=rwr_alpha,
            seed_quantile=seed_quantile
        )
        edge_weights_df = edge_weights_df.round(5)
        edge_weights_df.to_csv(os.path.join(save_path, f"{sample_id}_edges.csv"), index=False)
