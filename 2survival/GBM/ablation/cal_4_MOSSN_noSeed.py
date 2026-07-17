import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_uniform
import MOSSN_noSeed

available_datasets = ["CGGA301", "CGGA325", "CGGA693",
                      "GSE4412", "GSE13041", "GSE16011",
                      "GSE72951", "GSE74187", "TCGAGBM"]

researchAim = "GBM"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/MOSSN_noSeed"

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

    gamma = 2.0
    rwr_alpha = 0.3

    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = MOSSN_noSeed.MOSSN_noSeed_single_sample(
            sample_id=sample_id,
            G=G,
            uniform_original_weights=uniform_original_weights,
            expression_data=expression_data,
            gamma=gamma,
            rwr_alpha=rwr_alpha
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")
