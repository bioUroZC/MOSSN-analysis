import os
import sys

import pandas as pd
from tqdm import tqdm


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_NodeRWR


available_datasets = ["GSE10927", "GSE19750", "GSE33371",
                      "GSE49278", "GSE76019", "GSE76021",
                      "TCGAACC"]

researchAim = "ACC"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/NodeRWR"
    os.makedirs(save_path, exist_ok=True)

    links = pd.read_csv(
        PROJ_ROOT + "/1NT/1data/string/links.csv",
        index_col=0
    )
    expression_data = pd.read_csv(
        f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv",
        index_col=0
    )

    G, real_original_weights, expression_data = (
        Baseline_NodeRWR.prepare_data_NodeRWR(
            links=links,
            expression_data=expression_data
        )
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    gamma = 2.0
    rwr_alpha = 0.3
    seed_quantile = 0.9

    node_scores = {}
    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        node_scores[sample_id] = Baseline_NodeRWR.NodeRWR_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            gamma=gamma,
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
