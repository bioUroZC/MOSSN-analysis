import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ['GSE194040', 'GSE20194',  'GSE20271', 'GSE241876', 'GSE28844',
                      'GSE32646',  'GSE41998',  'GSE50948', 'GSE66305']

researchAim = "Paclitaxel"

dataset_organ_map = {
    'GSE194040': "Breast",
    'GSE20194': "Breast",
    'GSE20271': "Breast",
    'GSE241876': "Breast",
    'GSE28844': "Breast",
    'GSE32646': "Breast",
    'GSE41998': "Breast",
    'GSE50948': "Breast",
    'GSE66305': "Breast"
}

NormalFile = "/proj/c.zihao/work1/0ref/normal/combined_expr_df.csv"  


for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/3drugs/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/SSN"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    link_file = f"/proj/c.zihao/work1/1NT/1data/string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=dataset_organ_map[dataset_name],
        save_path=save_path
    )
