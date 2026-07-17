import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["GSE15459", "GSE26253",
                      "GSE26899", "GSE26901", "GSE29272", "GSE57303", "GSE62254",
                      "GSE84437", "TCGASTAD"]

researchAim = "STAD"

dataset_organ_map = {
    "GSE13861": "Stomach",
    "GSE15459": "Stomach",
    "GSE26253": "Stomach",
    "GSE26899": "Stomach",
    "GSE26901": "Stomach",          
    "GSE29272": "Stomach",
    "GSE57303": "Stomach",
    "GSE62254": "Stomach",
    "GSE84437": "Stomach",
    "TCGASTAD": "Stomach"   
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
