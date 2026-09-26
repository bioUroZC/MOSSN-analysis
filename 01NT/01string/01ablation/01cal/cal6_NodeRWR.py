import os
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import NodeRWR_func as Baseline_NodeRWR

available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

base_dir = "/proj/c.zihao/work1/01NT/01string/01ablation"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/NodeRWR"
    os.makedirs(save_path, exist_ok=True)

    links = pd.read_csv(
        "/proj/c.zihao/work1/01data/06string/links.csv",
        index_col=0
    )
    expression_data = pd.read_csv(
        f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv",
        index_col=0
    )

    G, real_original_weights, expression_data = (
        Baseline_NodeRWR.NodeRWR_prepare_data(
            links=links,
            expression_data=expression_data
        )
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    lam = 2.0
    rwr_alpha = 0.3
    seed_quantile = 0.9

    node_scores = {}
    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        node_scores[sample_id] = Baseline_NodeRWR.NodeRWR_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            lam=lam,
            rwr_alpha=rwr_alpha,
            seed_quantile=seed_quantile
        )

    result = pd.DataFrame(node_scores)
    result.index.name = "Interaction"
    result.to_csv(f"{save_path}/{dataset_name}.csv", float_format="%.5f")
    print(
        f"NodeRWR {dataset_name} -> "
        f"{result.shape[0]} nodes x {result.shape[1]} samples"
    )
