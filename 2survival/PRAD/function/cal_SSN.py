import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["DKFZ2018", "GSE116918", "GSE21034",
                      "GSE46602", "GSE54460", "GSE70768",
                      "GSE70769", "TCGAPRAD"]

researchAim = "PRAD"

dataset_organ_map = {
    "DKFZ2018": "Prostate",
    "GSE116918": "Prostate",
    "GSE21034": "Prostate",
    "GSE46602": "Prostate",
    "GSE54460": "Prostate",
    "GSE70768": "Prostate",
    "GSE70769": "Prostate",
    "TCGAPRAD": "Prostate",
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
