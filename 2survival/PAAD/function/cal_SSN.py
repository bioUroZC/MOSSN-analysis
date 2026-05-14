import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["CPTAC", "GSE28735", "GSE62452",
                      "GSE71729", "GSE79668", "GSE85916",
                      "MTAB6134", "QCMG", "TCGAPAAD"]

researchAim = "PAAD"

dataset_organ_map = {
    "CPTAC": "Pancreas",
    "GSE28735": "Pancreas",
    "GSE62452": "Pancreas",
    "GSE71729": "Pancreas",
    "GSE79668": "Pancreas",
    "GSE85916": "Pancreas",
    "MTAB6134": "Pancreas",
    "QCMG": "Pancreas",
    "TCGAPAAD": "Pancreas",
}

NormalFile = "/proj/c.zihao/work1/0ref/normal/combined_expr_df.csv"  


for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/SSN"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
