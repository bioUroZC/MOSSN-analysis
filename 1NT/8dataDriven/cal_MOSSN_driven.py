import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys
sys.path.append(r"//proj/c.zihao/work1/1NT/8dataDriven/")
import MOSSN_driven

researchAim = '8dataDriven'
cor_threshold = 0.9

available_datasets = ["BRCA", "CRC", "LIHC", "LUAD"]

base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/MOSSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    expression_data = pd.read_csv(
        f"/proj/c.zihao/work1/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv",
        index_col=0
    )

    corr_matrix, genes, expression_data = MOSSN_driven.prepare_data_MOSSN(
        expression_data=expression_data
    )

    G, real_original_weights, expression_data = MOSSN_driven.build_graph_from_corr(
        corr_matrix=corr_matrix,
        genes=genes,
        expression_data=expression_data,
        cor_threshold=cor_threshold
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    links_df = pd.DataFrame(
        [(u, v, w) for (u, v), w in real_original_weights.items()],
        columns=["protein1", "protein2", "score"]
    )
    links_df = links_df[["score", "protein1", "protein2"]]
    links_df.to_csv(f"{save_path}/links.csv", index=True)
    print(f"Links saved: {len(links_df)} edges")

    os.chdir(save_path)
    for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
        edge_weights_df = MOSSN_driven.MOSSN_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data
        )
        edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
