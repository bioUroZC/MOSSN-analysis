import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import PPIXpress

available_datasets = ["GSE102073", "GSE13876", "GSE140082",
                      "GSE17260", "GSE18520", "GSE23554",
                      "GSE26193", "GSE26712", "GSE30161",
                      "GSE31245", "GSE32062", "GSE51088",
                      "GSE53963", "GSE63885", "GSE73614",
                      "GSE8842", "GSE9891", "MTAB386",
                      "TCGAOV"]

researchAim = "OV"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/PPIXpress"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    PPIXpress.PPIXpressCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
