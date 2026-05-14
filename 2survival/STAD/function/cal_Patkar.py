import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import Patkar

available_datasets = ["GSE13861", "GSE15459", "GSE26253",
                      "GSE26899", "GSE26901", "GSE29272", "GSE57303", "GSE62254",
                      "GSE84437", "TCGASTAD"]

researchAim = "STAD"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/Patkar/"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    Patkar.PatkarCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
