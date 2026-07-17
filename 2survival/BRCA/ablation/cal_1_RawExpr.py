import os
import sys

import pandas as pd


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_RawExpr


available_datasets = ["GSE11121", "GSE12093", "GSE162228",
                      "GSE17705", "GSE20685", "GSE20711",
                      "GSE21653", "GSE22219", "GSE25055",
                      "GSE25065", "GSE42568", "GSE45255",
                      "GSE48390", "GSE61304", "GSE7390",
                      "TCGABRCA"]

researchAim = "BRCA"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/2survival/{researchAim}/"
    save_path = f"{base_dir}/{dataset_name}/RawExpr"
    os.makedirs(save_path, exist_ok=True)

    expression_data = pd.read_csv(
        f"{base_dir}/{dataset_name}/data/exprSet_filtered.csv",
        index_col=0
    )

    features = Baseline_RawExpr.get_raw_expression_features(expression_data)
    features.index.name = "Interaction"
    features.to_csv(f"{save_path}/{dataset_name}.csv", float_format="%.5f")
    print(
        f"RawExpr {dataset_name} -> "
        f"{features.shape[0]} genes x {features.shape[1]} samples"
    )
