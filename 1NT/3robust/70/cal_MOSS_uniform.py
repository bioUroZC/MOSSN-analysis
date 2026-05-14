import os
import pandas as pd
from tqdm import tqdm
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import MOSS_uniform

dataset_name = "LUAD"
researchAim = "3robust/70"
base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}"
save_path = f"{base_dir}/MOSS_uniform/{dataset_name}"
expr_file = f"{base_dir}/data/LUAD_exprSet_half.csv"

print(f"#========== {dataset_name} ==========")

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)
os.chdir(save_path)

seed_nodes_df = pd.read_csv("/proj/c.zihao/work1/1NT/1data/Markers/SSmarkers.csv", index_col=0)
links = pd.read_csv("/proj/c.zihao/work1/1NT/1data/string/links.csv", index_col=0)
expression_data = pd.read_csv(expr_file, index_col=0)

G, uniform_original_weights, string_original_weights, expression_data, seed_nodes_df = MOSS_uniform.prepare_data_MOSS_uniform(
    seed_nodes_df=seed_nodes_df,
    links=links,
    expression_data=expression_data
)
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

for sample_id in tqdm(expression_data.columns, desc="Processing Samples", unit="sample"):
    edge_weights_df = MOSS_uniform.MOSS_uniform_single_sample(
        sample_id=sample_id,
        G=G,
        uniform_original_weights=uniform_original_weights,
        string_original_weights=string_original_weights,
        expression_data=expression_data,
        seed_nodes_df=seed_nodes_df
    )
    edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
