import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Patkar

available_datasets = ["CGGA301", "CGGA325", "CGGA693",
                      "GSE4412", "GSE13041", "GSE16011",
                      "GSE72951", "GSE74187", "TCGAGBM"]

researchAim = "GBM"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/Patkar/"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/string/links.csv"

    Patkar.PatkarCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
