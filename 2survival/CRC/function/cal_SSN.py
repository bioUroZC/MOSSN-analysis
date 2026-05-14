import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

available_datasets = ["ACICAM", "GSE12945", "GSE17536",
                      "GSE17537", "GSE28722", "GSE29621",
                      "GSE39582", "GSE41258", "TCGACRC"]

researchAim = "CRC"

dataset_organ_map = {
    "ACICAM": "Colon",
    "GSE12945": "Colon",
    "GSE17536": "Colon",
    "GSE17537": "Colon",
    "GSE28722": "Colon",
    "GSE29621": "Colon",
    "GSE39582": "Colon",
    "GSE41258": "Colon",
    "TCGACRC": "Colon",
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
