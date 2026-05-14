"""MOSSN_Direct noise: direct cross-layer model with noisy EXP, LUAD only."""
import os, sys, shutil
import pandas as pd
from tqdm import tqdm

sys.path.append("/proj/c.zihao/work1/function/")
import MOSSN_coupled

CANCER       = "LUAD"
DATA_DIR     = "/proj/c.zihao/work1/4coupled/noise/noise_files"
LINK_PATH    = "/proj/c.zihao/work1/1NT/1data/string/links.csv"
OUT_ROOT     = "/proj/c.zihao/work1/4coupled/noise/results/MOSSN_Direct"
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
DIRECT_OMICS = ["MET", "CNV"]

links = pd.read_csv(LINK_PATH, index_col=0)

for k in NOISE_LEVELS:
    print(f"\n{'='*50}\nMOSSN_Direct  noise k={k}\n{'='*50}")
    save_dir = f"{OUT_ROOT}/k{k}/{CANCER}"
    if os.path.exists(save_dir): shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    omic_data = {
        "EXP": pd.read_csv(f"{DATA_DIR}/{CANCER}_EXP_k{k}.csv", index_col=0),
        "MET": pd.read_csv(f"{DATA_DIR}/{CANCER}_MET.csv", index_col=0),
        "CNV": pd.read_csv(f"{DATA_DIR}/{CANCER}_CNV.csv", index_col=0),
    }
    G, weights, omic_data, exp_genes = MOSSN_coupled.prepare_data_MOSSN_Direct(
        links, omic_data, direct_omics=DIRECT_OMICS)
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"EXP genes: {len(exp_genes)}, Samples: {omic_data['EXP'].shape[1]}")

    for sid in tqdm(omic_data["EXP"].columns, desc=f"k={k}"):
        df = MOSSN_coupled.MOSSN_Direct_single_sample(
            sid, G, weights, omic_data, exp_genes,
            direct_omics=DIRECT_OMICS, dynamic_cross=True)
        df["Node1"] = "EXP:" + df["Node1"]
        df["Node2"] = "EXP:" + df["Node2"]
        df.to_csv(f"{save_dir}/{sid}_edges.csv", index=False)
