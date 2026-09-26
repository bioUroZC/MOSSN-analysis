import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import SSN

dataset_name = "LUAD"
base_dir = "/proj/c.zihao/work1/02analysis/02noise"
LEVELS = ["05", "10", "15", "20"]
link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

for level in LEVELS:
    print(f"#========== {dataset_name} {level}% ==========")
    level_dir = f"{base_dir}/{level}"
    save_path = f"{level_dir}/SSN/{dataset_name}"
    NormalFile = f"{level_dir}/data/LUAD_normal_reference.csv"
    exprSetFile = f"{level_dir}/data/LUAD_exprSet_noise.csv"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=None,
        save_path=save_path
    )
