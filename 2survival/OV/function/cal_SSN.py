import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SSN

available_datasets = ["GSE102073", "GSE13876", "GSE17260", "GSE26193", "GSE26712", "GSE30161",
                      "GSE31245", "GSE32062", "GSE51088",
                      "GSE53963", "GSE8842", "GSE9891", "MTAB386",
                      "TCGAOV"]

researchAim = "OV"

dataset_organ_map = {
    "GSE102073": "Ovary",
    "GSE13876": "Ovary",
    "GSE140082": "Ovary",
    "GSE17260": "Ovary",
    "GSE18520": "Ovary",
    "GSE23554": "Ovary",
    "GSE26193": "Ovary",
    "GSE26712": "Ovary",
    "GSE30161": "Ovary",
    "GSE31245": "Ovary",
    "GSE32062": "Ovary",
    "GSE51088": "Ovary",
    "GSE53963": "Ovary",
    "GSE63885": "Ovary",
    "GSE73614": "Ovary",
    "GSE8842": "Ovary",
    "GSE9891": "Ovary",
    "MTAB386": "Ovary",
    "TCGAOV": "Ovary",
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
