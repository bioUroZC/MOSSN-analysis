import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import MOSSN_noPrior

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

researchAim = '2HuRI'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/"
    save_path = f"{base_dir}/MOSSN_noPrior/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    links = pd.read_csv(f"/proj/c.zihao/work1/1NT/1data/HuRI/links.csv", index_col=0)
    expression_data = pd.read_csv(f"/proj/c.zihao/work1/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv", index_col=0)

    G, real_original_weights, expression_data = MOSSN_noPrior.prepare_data_MOSSN_noPrior(
        links=links,
        expression_data=expression_data
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    gamma = 2.0
    rwr_alpha = 0.3
    seed_quantile = 0.9

    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = MOSSN_noPrior.MOSSN_noPrior_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            gamma=gamma,
            rwr_alpha=rwr_alpha,
            seed_quantile=seed_quantile
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
