import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import MOSSN_func as MOSSN
import NoRWR_func as EdgeNoRWR

available_datasets = ["CPTAC", "GSE167573", "TCGAKIRC"]

researchAim = "KIRC"

base_dir = f"/proj/c.zihao/work1/03survival/{researchAim}"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/{dataset_name}/noRWR"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(
        "/proj/c.zihao/work1/01data/06string/links.csv",
        index_col=0
    )
    expression_data = pd.read_csv(
        f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv",
        index_col=0
    )

    G, uniform_original_weights, expression_data = (
        MOSSN.MOSSN_prepare_data(
            links=links,
            expression_data=expression_data
        )
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    lam = 2.0

    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = EdgeNoRWR.NoRWR_single_sample(
            sample_id=sample_id,
            G=G,
            uniform_original_weights=uniform_original_weights,
            expression_data=expression_data,
            lam=lam
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")
