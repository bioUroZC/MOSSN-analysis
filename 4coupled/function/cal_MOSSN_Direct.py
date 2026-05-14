"""
MOSSN_Direct: primary model with direct same-gene cross-layer modulation.

Graph structure:
  EXP layer : TP53 ── MDM2 ── CDKN1A ...  (full PPI)
  MET layer :  MET:TP53    MET:MDM2  ...  (isolated nodes, no intra-layer edges)
  CNV layer :  CNV:TP53    CNV:MDM2  ...  (isolated nodes, no intra-layer edges)
  Cross-edges: TP53 <-> MET:TP53,  TP53 <-> CNV:TP53  (dynamic concordance weight)

MET/CNV values act as direct node-level modulators of EXP gene importance via RWR.
"""
import os, sys, shutil
import pandas as pd
from tqdm import tqdm

sys.path.append("/proj/c.zihao/work1/function/")
import MOSSN_coupled

CANCERS   = ["BLCA", "LIHC", "LUAD", "SARC", "STAD"]
DATA_DIR  = "/proj/c.zihao/work1/4coupled/files"

LINK_PATH = "/proj/c.zihao/work1/1NT/1data/string/links.csv"

OUT_ROOT  = "/proj/c.zihao/work1/4coupled/results/MOSSN_Direct"

DIRECT_OMICS = ["MET", "CNV"]

links = pd.read_csv(LINK_PATH, index_col=0)

for cancer in CANCERS:
    print(f"\n{'='*50}\n{cancer}\n{'='*50}")
    save_dir = f"{OUT_ROOT}/{cancer}"
    if os.path.exists(save_dir): shutil.rmtree(save_dir)
    os.makedirs(save_dir)

    omic_data = {
        "EXP": pd.read_csv(f"{DATA_DIR}/{cancer}_EXP.csv", index_col=0),
        "MET": pd.read_csv(f"{DATA_DIR}/{cancer}_MET.csv", index_col=0),
        "CNV": pd.read_csv(f"{DATA_DIR}/{cancer}_CNV.csv", index_col=0),
    }

    G, weights, omic_data, exp_genes = MOSSN_coupled.prepare_data_MOSSN_Direct(
        links, omic_data, direct_omics=DIRECT_OMICS)
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
          f"EXP genes: {len(exp_genes)}, Samples: {omic_data['EXP'].shape[1]}")

    for sid in tqdm(omic_data["EXP"].columns, desc=cancer):
        df = MOSSN_coupled.MOSSN_Direct_single_sample(
            sid, G, weights, omic_data, exp_genes,
            direct_omics=DIRECT_OMICS, dynamic_cross=True)
        df["Node1"] = "EXP:" + df["Node1"]
        df["Node2"] = "EXP:" + df["Node2"]
        df.to_csv(f"{save_dir}/{sid}_edges.csv", index=False)
