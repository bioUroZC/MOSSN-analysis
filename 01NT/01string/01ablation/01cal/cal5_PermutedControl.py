import os
import shutil
import sys

import pandas as pd

sys.path.append(r"/proj/c.zihao/work1/function/ablation")
import PermutedControl_func as Baseline_PermutedControl

available_datasets = [
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
]

base_dir = "/proj/c.zihao/work1/01NT/01string/01ablation"
benchmark_dir = "/proj/c.zihao/work1/01NT/01string/02benchmark/MOSSN"

for dataset_name in available_datasets:
    print(f"#========== {dataset_name} ==========")
    src_dir = f"{benchmark_dir}/{dataset_name}"
    save_path = f"{base_dir}/PermutedControl/{dataset_name}"

    if not os.path.isdir(src_dir):
        raise FileNotFoundError(
            f"Required benchmark result directory not found: {src_dir}. "
            "Run benchmark/MOSSN before PermutedControl."
        )

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)
    os.chdir(save_path)

    files = sorted(f for f in os.listdir(src_dir) if f.endswith("_edges.csv"))
    if not files:
        raise FileNotFoundError(
            f"No *_edges.csv files found in {src_dir}. "
            "PermutedControl requires completed benchmark/MOSSN outputs."
        )

    for f in files:
        sample_id = f[:-len("_edges.csv")]
        edge_df = pd.read_csv(os.path.join(src_dir, f))
        seed = Baseline_PermutedControl.seed_from_sample_id(sample_id, base_seed=1)
        permuted = Baseline_PermutedControl.permute_single_sample_edges(edge_df, seed=seed)
        permuted.to_csv(f"{sample_id}_edges.csv", index=False, float_format="%.5f")

    print(f"PermutedControl {dataset_name}: permuted {len(files)} sample edge files")
