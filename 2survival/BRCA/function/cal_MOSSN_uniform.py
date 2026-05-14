import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import MOSSN_uniform

available_datasets = ["GSE11121", "GSE12093", "GSE162228",
                      "GSE17705", "GSE20685", "GSE20711",
                      "GSE21653", "GSE22219", "GSE25055",
                      "GSE25065", "GSE42568", "GSE45255",
                      "GSE48390", "GSE61304", "GSE7390",
                      "TCGABRCA"]

researchAim = "BRCA"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/MOSSN_uniform"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(f"/proj/c.zihao/work1/1NT/1data/string/links.csv", index_col=0)
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
