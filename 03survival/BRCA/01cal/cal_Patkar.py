import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import Patkar

available_datasets = [
    "GSE11121", "GSE162228", "GSE17705", "GSE20685",
    "GSE20711", "GSE21653", "GSE22219", "GSE25055", "GSE25065",
    "GSE45255", "GSE61304", "GSE7390", "TCGABRCA",
]

research_aim = "BRCA"
base_dir = f"/proj/c.zihao/work1/03survival/{research_aim}"
link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/{dataset_name}/Patkar"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv"
    if not os.path.isfile(exprSetFile):
        raise FileNotFoundError(
            f"Filtered expression matrix not found for {dataset_name}: "
            f"{exprSetFile}"
        )

    Patkar.PatkarCal(
        exprSetFile=exprSetFile,
        link_file=link_file,
        save_path=save_path
    )
