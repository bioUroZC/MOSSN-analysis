import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["GSE30219", "GSE31210", "GSE41271",
                      "GSE42127", "GSE50081", "GSE68465",
                      "GSE72094", "TCGALUAD"]

researchAim = "LUAD"

dataset_organ_map = {
    "GSE30219": "Lung",
    "GSE31210": "Lung",
    "GSE41271": "Lung",
    "GSE42127": "Lung",
    "GSE50081": "Lung",
    "GSE68465": "Lung",
    "GSE72094": "Lung",
    "TCGALUAD": "Lung",
}

NormalFile = "/proj/c.zihao/work1/0ref/normal/combined_expr_df.csv"  


for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"/proj/c.zihao/work1/2survival/{researchAim}/"
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
