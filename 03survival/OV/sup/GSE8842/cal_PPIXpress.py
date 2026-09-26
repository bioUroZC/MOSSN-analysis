import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import PPIXpress

available_datasets = ["GSE8842"]

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
