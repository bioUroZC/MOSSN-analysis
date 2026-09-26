import argparse
import os
import shutil
import sys

import pandas as pd
sys.path.append(r"/proj/c.zihao/work1/function/benchmark/")
import SSN

available_datasets = ["LUAD", "eLUAD", "KIRC", "PRAD"]
dataset_dirs = {"LUAD": "01LUAD", "eLUAD": "02eLUAD", "KIRC": "03KIRC", "PRAD": "04PRAD"}
normal_files = {
    "LUAD": "/proj/c.zihao/work1/01NT/03GEO/01LUAD/LUAD_normal_exprSet_filtered.csv",
    "eLUAD": "/proj/c.zihao/work1/01NT/03GEO/02eLUAD/eLUAD_normal_exprSet_filtered.csv",
    "KIRC": "/proj/c.zihao/work1/01NT/03GEO/03KIRC/KIRC_normal_exprSet_filtered.csv",
    "PRAD": "/proj/c.zihao/work1/01NT/03GEO/04PRAD/PRAD_normal_exprSet_filtered.csv",
}

parser = argparse.ArgumentParser(
    description="Calculate SSN networks using each GEO cohort's own normal samples."
)
parser.add_argument(
    "--dataset",
    choices=available_datasets,
    help="Run only one dataset. Omit this option to run all datasets serially."
)
args = parser.parse_args()
datasets_to_run = [args.dataset] if args.dataset else available_datasets

base_dir = "/proj/c.zihao/work1/01NT/03GEO/Output"

for dataset_name in datasets_to_run:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/SSN/{dataset_name}"

    exprSetFile = f"/proj/c.zihao/work1/01NT/03GEO/{dataset_dirs[dataset_name]}/{dataset_name}_exprSet_filtered.csv"
    NormalFile = normal_files[dataset_name]
    if not os.path.isfile(NormalFile):
        raise FileNotFoundError(
            f"Cohort-matched normal reference not found: {NormalFile}. Run "
            f"/proj/c.zihao/work1/01NT/03GEO/{dataset_dirs[dataset_name]}/2normal.R first."
        )

    normal_reference = pd.read_csv(NormalFile, index_col=0)
    if normal_reference.shape[1] < 3:
        raise ValueError(
            f"{dataset_name} has only {normal_reference.shape[1]} normal reference samples; "
            "at least 3 are required because normal samples are evaluated "
            "with leave-one-out SSN references."
        )
    expression_columns = pd.read_csv(exprSetFile, index_col=0, nrows=0).columns
    missing_samples = normal_reference.columns.difference(expression_columns)
    if len(missing_samples) > 0:
        raise ValueError(
            f"Normal reference samples absent from {exprSetFile}: "
            f"{', '.join(missing_samples)}"
        )

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    print(
        f"Using {normal_reference.shape[1]} cohort-matched normal samples as the "
        f"SSN reference: {NormalFile}"
    )

    link_file = "/proj/c.zihao/work1/01data/06string/links.csv"

    SSN.SSNcal(
        exprSetFile=exprSetFile,
        NormalFile=NormalFile,
        link_file=link_file,
        organ=None,
        save_path=save_path
    )
