import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_uniform
import pandas as pd
from tqdm import tqdm

dataset_name = "LUAD"
researchAim = "5noise/10"
base_dir = f"{PROJ_ROOT}/1NT/{researchAim}"
save_path = f"{base_dir}/MOSSN_uniform/{dataset_name}"

link_file = PROJ_ROOT + "/1NT/1data/string/links.csv"
expr_file = f"{base_dir}/data/LUAD_exprSet_noise.csv"

print(f"#========== {dataset_name} ==========")

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)
os.chdir(save_path)

links = pd.read_csv(link_file, index_col=0)
expression_data = pd.read_csv(expr_file, index_col=0)

G, uniform_original_weights, expression_data = MOSSN_uniform.prepare_data_MOSSN_uniform(
    links=links,
    expression_data=expression_data
)
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
    edge_weights_df.to_csv(f"{sample_id}_edges.csv", index=False)
