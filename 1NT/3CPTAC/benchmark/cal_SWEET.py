import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import SWEET

available_datasets = ["LUAD", "eLUAD", "KIRC", "PRAD"]

researchAim = '2CPTAC'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{researchAim}/"
    save_path = f"{base_dir}/SWEET/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{PROJ_ROOT}/1NT/2CPTAC/{dataset_name}/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/string/links.csv"

    SWEET.SWEETcal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
