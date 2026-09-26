import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import PPIXpress

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

base_dir = "/proj/c.zihao/work1/02analysis/03net/2intact/out"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/PPIXpress/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = "/proj/c.zihao/work1/01data/07intact/intact_link.csv"

    PPIXpress.PPIXpressCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
