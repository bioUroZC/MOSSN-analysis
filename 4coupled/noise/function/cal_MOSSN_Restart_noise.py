"""MOSSN_Restart noise: restart-coupled variant with noisy EXP, LUAD only."""
import os, sys, shutil
import pandas as pd
from tqdm import tqdm


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_coupled

CANCER       = "LUAD"
DATA_DIR     = PROJ_ROOT + "/4coupled/noise/noise_files"
LINK_PATH    = PROJ_ROOT + "/1NT/1data/string/links.csv"
OUT_ROOT     = PROJ_ROOT + "/4coupled/noise/results/MOSSN_Restart"
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

links = pd.read_csv(LINK_PATH, index_col=0)

for k in NOISE_LEVELS:
    print(f"\n{'='*50}\nMOSSN_Restart  noise k={k}\n{'='*50}")
    save_dir = f"{OUT_ROOT}/k{k}/{CANCER}"
    if os.path.exists(save_dir): shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    omic_data = {
        "EXP": pd.read_csv(f"{DATA_DIR}/{CANCER}_EXP_k{k}.csv", index_col=0),
        "MET": pd.read_csv(f"{DATA_DIR}/{CANCER}_MET.csv", index_col=0),
        "CNV": pd.read_csv(f"{DATA_DIR}/{CANCER}_CNV.csv", index_col=0),
    }
    G, weights, omic_data = MOSSN_coupled.prepare_data_MOSSN_Coupled(links, omic_data)
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"Samples: {omic_data['EXP'].shape[1]}")

    for sid in tqdm(omic_data["EXP"].columns, desc=f"k={k}"):
        df = MOSSN_coupled.MOSSN_Coupled_single_sample(
            sid, G, weights, omic_data, coupled_omics=["MET", "CNV"],
            dynamic_cross=True, beta=0.7)
        df["Node1"] = "EXP:" + df["Node1"]
        df["Node2"] = "EXP:" + df["Node2"]
        df.to_csv(f"{save_dir}/{sid}_edges.csv", index=False)
