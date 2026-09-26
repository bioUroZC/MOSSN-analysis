import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import Patkar

available_datasets = ["LUAD", "eLUAD", "KIRC", "PRAD"]
dataset_dirs = {"LUAD": "01LUAD", "eLUAD": "02eLUAD", "KIRC": "03KIRC", "PRAD": "04PRAD"}

base_dir = "/proj/c.zihao/work1/01NT/03GEO/Output"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/Patkar/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/01NT/03GEO/{dataset_dirs[dataset_name]}/{dataset_name}_exprSet_filtered.csv"
    link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

    Patkar.PatkarCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
