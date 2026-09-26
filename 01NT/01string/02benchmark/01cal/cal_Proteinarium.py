import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import Proteinarium

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

base_dir = "/proj/c.zihao/work1/01NT/01string/02benchmark"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/Proteinarium/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

    Proteinarium.ProteinariumCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
