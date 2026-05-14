import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SWEET

available_datasets = ['GSE194040', 'GSE20194', 'GSE20271', 'GSE241876', 'GSE28844',
                      'GSE32646', 'GSE41998', 'GSE50948', 'GSE66305']

researchAim = "Paclitaxel"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/3drugs/{researchAim}/"
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
