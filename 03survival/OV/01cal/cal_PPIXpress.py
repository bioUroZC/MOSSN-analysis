import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import PPIXpress

available_datasets = ["GSE102073", "GSE26193", "GSE26712", "GSE30161",
                      "GSE31245", "GSE32062", "GSE51088", "GSE53963",
                      "GSE8842", "GSE9891", "MTAB386", "TCGAOV"]

researchAim = "OV"

base_dir = f"/proj/c.zihao/work1/03survival/{researchAim}"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/{dataset_name}/PPIXpress"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

    PPIXpress.PPIXpressCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
