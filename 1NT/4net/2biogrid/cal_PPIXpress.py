import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import PPIXpress

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

researchAim = '4net/2biogrid'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{researchAim}/"
    save_path = f"{base_dir}/PPIXpress/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{PROJ_ROOT}/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/biogrid/biogrid_link.csv"

    PPIXpress.PPIXpressCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
