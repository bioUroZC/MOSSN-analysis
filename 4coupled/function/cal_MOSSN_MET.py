"""MOSSN_MET: single-layer RWR on methylation (using the EXP slot)."""
import os, sys, shutil
import pandas as pd
from tqdm import tqdm


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_coupled

CANCERS   = ["ACC", "BLCA", "BRCA", "CESC", 
             "CRC", "ESCA", "GBM", "HNSC", "KIRC",
             "LGG", "LIHC", "LUAD", "LUSC", 
             "PAAD", "PRAD", "SARC", "STAD"]
DATA_DIR  = PROJ_ROOT + "/4coupled/files"

LINK_PATH = PROJ_ROOT + "/1NT/1data/string/links.csv"

OUT_ROOT  = PROJ_ROOT + "/4coupled/results/MOSSN_MET"

links = pd.read_csv(LINK_PATH, index_col=0)

for cancer in CANCERS:
    print(f"\n{'='*50}\n{cancer}\n{'='*50}")
    save_dir = f"{OUT_ROOT}/{cancer}"
    if os.path.exists(save_dir): shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    # use MET as the "EXP" slot so single-layer RWR runs on methylation
    omic_data = {"EXP": pd.read_csv(f"{DATA_DIR}/{cancer}_MET.csv", index_col=0)}
    G, weights, omic_data = MOSSN_coupled.prepare_data_MOSSN_Coupled(links, omic_data)
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, Samples: {omic_data['EXP'].shape[1]}")

    for sid in tqdm(omic_data["EXP"].columns, desc=cancer):
        df = MOSSN_coupled.MOSSN_Coupled_single_sample(
            sid, G, weights, omic_data, coupled_omics=[], beta=1.0)
        df["Node1"] = "MET:" + df["Node1"]
        df["Node2"] = "MET:" + df["Node2"]
        df.to_csv(f"{save_dir}/{sid}_edges.csv", index=False)
