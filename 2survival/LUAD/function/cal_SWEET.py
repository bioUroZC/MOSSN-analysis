import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SWEET

available_datasets = ["GSE30219", "GSE31210", "GSE41271",
                      "GSE42127", "GSE50081", "GSE68465",
                      "GSE72094", "TCGALUAD"]

researchAim = "LUAD"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/SWEET"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    SWEET.SWEETcal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
