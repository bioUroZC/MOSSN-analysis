import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import PPIXpress

available_datasets = ["DKFZ2018", "GSE116918", "GSE21034",
                      "GSE46602", "GSE54460", "GSE70768",
                      "GSE70769", "TCGAPRAD"]

researchAim = "PRAD"

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
