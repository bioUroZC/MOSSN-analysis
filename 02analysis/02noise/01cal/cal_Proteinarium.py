import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import Proteinarium

dataset_name = "LUAD"
base_dir = "/proj/c.zihao/work1/02analysis/02noise"
LEVELS = ["05", "10", "15", "20"]
link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

for level in LEVELS:
    print(f"#========== {dataset_name} {level}% ==========")
    level_dir = f"{base_dir}/{level}"
    save_path = f"{level_dir}/Proteinarium/{dataset_name}"
    exprSetFile = f"{level_dir}/data/LUAD_exprSet_noise.csv"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    Proteinarium.ProteinariumCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
