"""
MOSSN_MultiLayer: full multilayer comparative variant on EXP/MET/CNV.

Graph structure:
  - EXP, MET, and CNV each keep their own full PPI layer
  - The same gene is connected across all omic layers
  - One random walk is performed on the whole multilayer graph

This compares the primary direct model against a full multilayer propagation
scheme where every omic keeps its own PPI layer.
"""
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

OUT_ROOT  = PROJ_ROOT + "/4coupled/results/MOSSN_MultiLayer"

ACTIVE_OMICS = ["EXP", "MET", "CNV"]

links = pd.read_csv(LINK_PATH, index_col=0)

for cancer in CANCERS:
    print(f"\n{'='*50}\n{cancer}\n{'='*50}")
    save_dir = f"{OUT_ROOT}/{cancer}"
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    omic_data = {
        "EXP": pd.read_csv(f"{DATA_DIR}/{cancer}_EXP.csv", index_col=0),
        "MET": pd.read_csv(f"{DATA_DIR}/{cancer}_MET.csv", index_col=0),
        "CNV": pd.read_csv(f"{DATA_DIR}/{cancer}_CNV.csv", index_col=0),
    }

    G, weights, omic_data, exp_genes = MOSSN_coupled.prepare_data_MOSSN_Multilayer(
        links, omic_data, active_omics=ACTIVE_OMICS
    )
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"EXP genes: {len(exp_genes)}, Samples: {omic_data['EXP'].shape[1]}")

    for sid in tqdm(omic_data["EXP"].columns, desc=cancer):
        df = MOSSN_coupled.MOSSN_Multilayer_single_sample(
            sid, G, weights, omic_data, exp_genes,
            active_omics=ACTIVE_OMICS, dynamic_cross=True
        )
        df["Node1"] = "EXP:" + df["Node1"]
        df["Node2"] = "EXP:" + df["Node2"]
        df.to_csv(f"{save_dir}/{sid}_edges.csv", index=False)
