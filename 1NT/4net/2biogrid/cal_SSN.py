import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

dataset_organ_map = {
    "BLCA": "Bladder",
    "BRCA": "Breast",
    "CRC":  "Colon",
    "ESCA": "Esophagus",
    "HNSC": "Salivary Gland",  
    "KIRC": "Kidney",
    "LIHC": "Liver",
    "LUAD": "Lung",
    "LUSC": "Lung",
    "PRAD": "Prostate",
    "STAD": "Stomach",
}

NormalFile = PROJ_ROOT + "/0ref/GTEx/combined_expr_df.csv"  
researchAim = '4net/2biogrid'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{researchAim}/"
    save_path = f"{base_dir}/SSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{PROJ_ROOT}/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/biogrid/biogrid_link.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
