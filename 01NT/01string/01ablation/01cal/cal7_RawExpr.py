import os
import sys

import pandas as pd

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import RawExpr_func as Baseline_RawExpr

available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

base_dir = "/proj/c.zihao/work1/01NT/01string/01ablation"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    save_path = f"{base_dir}/RawExpr"
    os.makedirs(save_path, exist_ok=True)

    expression_data = pd.read_csv(
        f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv",
        index_col=0
    )

    features = Baseline_RawExpr.get_raw_expression_features(expression_data)
    features.index.name = "Interaction"
    features.to_csv(f"{save_path}/{dataset_name}.csv", float_format="%.5f")
    print(
        f"RawExpr {dataset_name} -> "
        f"{features.shape[0]} genes x {features.shape[1]} samples"
    )
