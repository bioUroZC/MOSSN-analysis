import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["GSE107943", "GSE244807", "GSE89748",
                      "MTAB6389", "TCGACHOL"]

researchAim = "CHOL"

dataset_organ_map = {
    "GSE107943": "Liver",
    "GSE244807": "Liver",
    "GSE89748": "Liver",
    "MTAB6389": "Liver",
    "TCGACHOL": "Liver",
}

NormalFile = PROJ_ROOT + "/0ref/GTEx/combined_expr_df.csv"  


for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/SSN"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
