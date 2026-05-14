import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import PPIXpress

available_datasets = ["CGGA301", "CGGA325", "CGGA693",
                      "GSE16011", "MTAB3892", "TCGALGG"]

researchAim = "LGG"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/PPIXpress"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    PPIXpress.PPIXpressCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
