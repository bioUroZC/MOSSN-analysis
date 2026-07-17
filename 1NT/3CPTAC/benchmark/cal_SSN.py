import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["LUAD", "eLUAD", "KIRC", "PRAD"]

researchAim = '2CPTAC'

dataset_organ_map = {
    "LUAD": "Lung",
    "eLUAD": "Lung",
    "KIRC": "Kidney",
    "PRAD": "Prostate"
}

NormalFile = PROJ_ROOT + "/0ref/GTEx/combined_expr_df.csv"  
researchAim = '2CPTAC'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{researchAim}/"
    save_path = f"{base_dir}/SSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{PROJ_ROOT}/1NT/2CPTAC/{dataset_name}/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
