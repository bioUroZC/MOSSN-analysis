import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import RandomBackbone_func as Baseline_RandomBackbone
import MOSSN_func as MOSSN

available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

base_dir = "/proj/c.zihao/work1/01NT/01string/01ablation"

for dataset_seed, dataset_name in enumerate(available_datasets, start=1):
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/RandomBackbone/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(
        "/proj/c.zihao/work1/01data/06string/links.csv",
        index_col=0
    )
    expression_data = pd.read_csv(
        f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv",
        index_col=0
    )

    G, _, expression_data = Baseline_RandomBackbone.RandomBackbone_prepare_data(
        links=links,
        expression_data=expression_data,
        seed=dataset_seed
    )
    uniform_original_weights = {(u, v): 1.0 for u, v in G.edges()}
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

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
            seed_quantile=seed_quantile
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")
