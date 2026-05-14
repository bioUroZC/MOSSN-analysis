import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import Proteinarium

available_datasets = ["GSE116174", "GSE144269", "GSE14520",
                      "GSE54236", "ICGC", "TCGALIHC"]

researchAim = "LIHC"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/Proteinarium"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    Proteinarium.ProteinariumCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
