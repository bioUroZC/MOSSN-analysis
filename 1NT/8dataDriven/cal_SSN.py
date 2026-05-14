import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["BRCA", "CRC", "LIHC", "LUAD"]

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

NormalFile = "/proj/c.zihao/work1/0ref/normal/combined_expr_df.csv"  
researchAim = '8dataDriven'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/"
    save_path = f"{base_dir}/SSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{base_dir}/MOSSN/{dataset_name}/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
