import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["CGGA301", "CGGA325", "CGGA693",
                      "GSE4412", "GSE13041", "GSE16011",
                      "GSE72951", "GSE74187", "TCGAGBM"]

researchAim = "GBM"

dataset_organ_map = {
    "CGGA301": "Brain",
    "CGGA325": "Brain",
    "CGGA693": "Brain",
    "GSE4412": "Brain",
    "GSE13041": "Brain",
    "GSE16011": "Brain",
    "GSE72951": "Brain",
    "GSE74187": "Brain",
    "GSE83300": "Brain",
    "TCGAGBM": "Brain",
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
