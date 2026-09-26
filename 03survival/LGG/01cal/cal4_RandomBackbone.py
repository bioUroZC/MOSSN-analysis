import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import RandomBackbone_func as Baseline_RandomBackbone
import MOSSN_func as MOSSN

available_datasets = ["CGGA301", "CGGA325",
                      "GSE16011", "TCGALGG"]

researchAim = "LGG"

base_dir = f"/proj/c.zihao/work1/03survival/{researchAim}"

SEED = 1

links = pd.read_csv(
    "/proj/c.zihao/work1/01data/06string/links.csv",
    index_col=0
)

G, uniform_original_weights = Baseline_RandomBackbone.RandomBackbone_prepare_shared(
    links=links,
    seed=SEED
)
print(f"Shared random backbone -- Nodes: {G.number_of_nodes()}, "
      f"Edges: {G.number_of_edges()}")

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/{dataset_name}/RandomBackbone"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    expression_data = pd.read_csv(
        f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv",
        index_col=0
    )
    expression_data = expression_data.loc[expression_data.index.isin(G.nodes)]

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
