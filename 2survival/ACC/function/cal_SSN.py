import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["GSE10927", "GSE19750", "GSE33371",
                      "GSE49278", "GSE76019", "GSE76021",
                      "TCGAACC"]

researchAim = "ACC"

dataset_organ_map = {
    "GSE10927": "Adrenal Gland",
    "GSE19750": "Adrenal Gland",
    "GSE33371": "Adrenal Gland",
    "GSE49278": "Adrenal Gland",
    "GSE76019": "Adrenal Gland",
    "GSE76021": "Adrenal Gland",
    "TCGAACC": "Adrenal Gland",
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
