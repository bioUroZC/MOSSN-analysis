import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_RandomBackbone
import MOSSN_uniform


available_datasets = ["GSE11121", "GSE12093", "GSE162228",
                      "GSE17705", "GSE20685", "GSE20711",
                      "GSE21653", "GSE22219", "GSE25055",
                      "GSE25065", "GSE42568", "GSE45255",
                      "GSE48390", "GSE61304", "GSE7390",
                      "TCGABRCA"]

researchAim = "BRCA"

for dataset_seed, dataset_name in enumerate(available_datasets, start=1):
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/RandomBackbone"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(
        PROJ_ROOT + "/1NT/1data/string/links.csv",
        index_col=0
    )
    expression_data = pd.read_csv(
        f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv",
        index_col=0
    )

    G, _, expression_data = Baseline_RandomBackbone.prepare_data_RandomBackbone(
        links=links,
        expression_data=expression_data,
        seed=dataset_seed
    )
    uniform_original_weights = {(u, v): 1.0 for u, v in G.edges()}
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    gamma = 2.0
    rwr_alpha = 0.3
    seed_quantile = 0.9

    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = MOSSN_uniform.MOSSN_uniform_single_sample(
            sample_id=sample_id,
            G=G,
            uniform_original_weights=uniform_original_weights,
            expression_data=expression_data,
            gamma=gamma,
            rwr_alpha=rwr_alpha,
            seed_quantile=seed_quantile
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")
