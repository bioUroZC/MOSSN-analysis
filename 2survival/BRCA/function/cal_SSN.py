import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["GSE11121", "GSE12093", "GSE162228",
                      "GSE17705", "GSE20685", "GSE20711",
                      "GSE21653", "GSE22219", "GSE25055",
                      "GSE25065", "GSE42568", "GSE45255",
                      "GSE48390", "GSE61304", "GSE7390",
                      "TCGABRCA"]

researchAim = "BRCA"

dataset_organ_map = {
    "GSE11121": "Breast",
    "GSE12093": "Breast",
    "GSE162228": "Breast",
    "GSE17705": "Breast",
    "GSE20685": "Breast",
    "GSE20711": "Breast",
    "GSE21653": "Breast",
    "GSE22219": "Breast",
    "GSE25055": "Breast",
    "GSE25065": "Breast",
    "GSE42568": "Breast",
    "GSE45255": "Breast",
    "GSE48390": "Breast",
    "GSE61304": "Breast",
    "GSE7390": "Breast",
    "TCGABRCA": "Breast",
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
