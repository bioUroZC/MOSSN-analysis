import os
import sys

import pandas as pd


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_RawExpr


available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

research_aim = "2string"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{research_aim}/ablation"
    save_path = f"{base_dir}/RawExpr"
    os.makedirs(save_path, exist_ok=True)

    expression_data = pd.read_csv(
        f"{PROJ_ROOT}/1NT/1data/exprset/{dataset_name}_exprSet_filtered.csv",
        index_col=0
    )

    features = Baseline_RawExpr.get_raw_expression_features(expression_data)
    features.index.name = "Interaction"
    features.to_csv(f"{save_path}/{dataset_name}.csv", float_format="%.5f")
    print(
        f"RawExpr {dataset_name} -> "
        f"{features.shape[0]} genes x {features.shape[1]} samples"
    )
