import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ['GSE140494', 'GSE20194', 'GSE20271', 'GSE22093', 'GSE23988',
                      'GSE32646', 'GSE42822', 'GSE50948', 'GSE66305', 'GSE6861']

researchAim = "FluoroBRCA"

dataset_organ_map = {
    'GSE140494': "Breast", 'GSE20194': "Breast", 'GSE20271': "Breast",
    'GSE22093':  "Breast", 'GSE23988': "Breast", 'GSE32646': "Breast",
    'GSE42822':  "Breast", 'GSE50948': "Breast", 'GSE66305': "Breast",
    'GSE6861':   "Breast"
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
