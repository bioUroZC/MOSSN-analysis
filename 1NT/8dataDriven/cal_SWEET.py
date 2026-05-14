import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SWEET

available_datasets = ["BRCA", "CRC", "LIHC", "LUAD"]

researchAim = '8dataDriven'

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/"
    save_path = f"{base_dir}/SWEET/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = f"{base_dir}/MOSSN/{dataset_name}/links.csv"

    SWEET.SWEETcal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
