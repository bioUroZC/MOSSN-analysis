import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["GSE13507", "GSE31684", "GSE32894",
                      "GSE48276", "TCGABLCA"]

researchAim = "BLCA"

dataset_organ_map = {
    "GSE13507": "Bladder",
    "GSE31684": "Bladder",
    "GSE32894": "Bladder",
    "GSE48276": "Bladder",
    "TCGABLCA": "Bladder",
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
