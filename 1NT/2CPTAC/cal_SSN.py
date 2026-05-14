import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["LUAD", "eLUAD", "KIRC", "PRAD"]

researchAim = '2CPTAC'

dataset_organ_map = {
    "LUAD": "Lung",
    "eLUAD": "Lung",
    "KIRC": "Kidney",
    "PRAD": "Prostate"
}

NormalFile = "/proj/c.zihao/work1/0ref/normal/combined_expr_df.csv"  
researchAim = '2CPTAC'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/"
    save_path = f"{base_dir}/SSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/1NT/2CPTAC/{dataset_name}/{dataset_name}_exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
