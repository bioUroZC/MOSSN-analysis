import argparse
import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import SSN

available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                      "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

parser = argparse.ArgumentParser(
    description="Calculate SSN networks using cancer-specific paired-normal references."
)
parser.add_argument(
    "--cancer",
    choices=available_datasets,
    help="Run only one cancer type. Omit this option to run all cancer types serially."
)
args = parser.parse_args()
datasets_to_run = [args.cancer] if args.cancer else available_datasets

base_dir = "/proj/c.zihao/work1/01NT/01string/02benchmark"
normal_dir = "/proj/c.zihao/work1/01data/03TCGAnormal/paired_normals"

for dataset_name in datasets_to_run:
    print(f"#========== {dataset_name} ==========")
    
    normal_file = f"{normal_dir}/TCGA-{dataset_name}_paired_normal.csv"
    if not os.path.isfile(normal_file):
        raise FileNotFoundError(
            f"Paired-normal reference not found: {normal_file}. Run "
            "/proj/c.zihao/work1/01data/03TCGAnormal/01_extract_paired_normals.R first."
        )

    save_path = f"{base_dir}/SSN/{dataset_name}"

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    exprSetFile = f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv"
    link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=normal_file,
        link_file=link_file,
        organ=None,
        save_path=save_path
    )
