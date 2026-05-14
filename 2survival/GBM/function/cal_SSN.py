import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["CGGA301", "CGGA325", "CGGA693",
                      "GSE4412", "GSE13041", "GSE16011",
                      "GSE72951", "GSE74187", "GSE83300",
                      "TCGAGBM"]

researchAim = "GBM"

dataset_organ_map = {
    "CGGA301": "Brain",
    "CGGA325": "Brain",
    "CGGA693": "Brain",
    "GSE4412": "Brain",
    "GSE13041": "Brain",
    "GSE16011": "Brain",
    "GSE72951": "Brain",
    "GSE74187": "Brain",
    "GSE83300": "Brain",
    "TCGAGBM": "Brain",
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
