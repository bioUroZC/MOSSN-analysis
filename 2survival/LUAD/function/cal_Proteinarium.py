import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Proteinarium

available_datasets = ["GSE30219", "GSE31210", "GSE41271",
                      "GSE42127", "GSE50081", "GSE68465",
                      "GSE72094", "TCGALUAD"]

researchAim = "LUAD"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/Proteinarium"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"{PROJ_ROOT}/1NT/1data/string/links.csv"

    Proteinarium.ProteinariumCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
