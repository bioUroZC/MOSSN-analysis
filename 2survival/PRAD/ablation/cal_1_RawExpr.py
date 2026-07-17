import os
import sys

import pandas as pd


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_RawExpr


available_datasets = ["DKFZ2018", "GSE116918", "GSE21034",
                      "GSE46602", "GSE54460", "GSE70768",
                      "GSE70769", "TCGAPRAD"]

researchAim = "PRAD"

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
