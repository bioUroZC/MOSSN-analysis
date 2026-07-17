import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_uniform

available_datasets = ["CPTAC", "GSE28735", "GSE62452",
                      "GSE71729", "GSE79668", "GSE85916",
                      "MTAB6134", "QCMG", "TCGAPAAD"]

researchAim = "PAAD"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/MOSSN_uniform"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(f"{PROJ_ROOT}/1NT/1data/string/links.csv", index_col=0)
    expression_data = pd.read_csv(f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv", index_col=0)

    G, uniform_original_weights, expression_data = MOSSN_uniform.prepare_data_MOSSN_uniform(
        links=links,
        expression_data=expression_data
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = MOSSN_uniform.MOSSN_uniform_single_sample(
            sample_id=sample_id,
            G=G,
            uniform_original_weights=uniform_original_weights,
            expression_data=expression_data
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
