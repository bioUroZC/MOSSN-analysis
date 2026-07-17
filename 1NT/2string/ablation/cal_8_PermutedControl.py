import os
import shutil
import sys

import pandas as pd


PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import Baseline_PermutedControl


available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

research_aim = "2string"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    base_dir = f"{PROJ_ROOT}/1NT/{research_aim}/ablation"
    src_dir = f"{PROJ_ROOT}/1NT/{research_aim}/benchmark/MOSSN_uniform/{dataset_name}"
    save_path = f"{base_dir}/PermutedControl/{dataset_name}"

    if not os.path.isdir(src_dir):
        raise FileNotFoundError(
            f"Required benchmark result directory not found: {src_dir}. "
            "Run benchmark/MOSSN_uniform before PermutedControl."
        )

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    files = sorted(f for f in os.listdir(src_dir) if f.endswith("_edges.csv"))
    if not files:
        raise FileNotFoundError(
            f"No *_edges.csv files found in {src_dir}. "
            "PermutedControl requires completed benchmark/MOSSN_uniform outputs."
        )

    for f in files:
        sample_id = f[:-len("_edges.csv")]
        edge_df = pd.read_csv(os.path.join(src_dir, f))
        seed = Baseline_PermutedControl.seed_from_sample_id(sample_id, base_seed=1)
        permuted = Baseline_PermutedControl.permute_single_sample_edges(edge_df, seed=seed)
        permuted.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")

    print(f"PermutedControl {dataset_name}: permuted {len(files)} sample edge files")
