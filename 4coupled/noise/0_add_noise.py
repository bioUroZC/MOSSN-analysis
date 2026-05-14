"""
Generate noisy EXP data for LUAD at multiple noise levels.

Only the EXP matrix is perturbed:
    noise_std = k * per-gene std(EXP across samples)

MET, CNV, and OS are copied unchanged into noise_files/ for convenience.
"""
import os, shutil
import numpy as np
import pandas as pd

CANCER      = "LUAD"
DATA_DIR    = "/proj/c.zihao/work1/4coupled/files"
OUT_DIR     = "/proj/c.zihao/work1/4coupled/noise/noise_files"
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
SEED        = 42

os.makedirs(OUT_DIR, exist_ok=True)

exp = pd.read_csv(f"{DATA_DIR}/{CANCER}_EXP.csv", index_col=0)
gene_stds = exp.std(axis=1)
rng = np.random.default_rng(SEED)

for k in NOISE_LEVELS:
    if k == 0.0:
        noisy = exp.copy()
    else:
        noise = rng.normal(0, 1, size=exp.shape) * (gene_stds.values[:, None] * k)
        noisy = exp + noise
    out_path = f"{OUT_DIR}/{CANCER}_EXP_k{k}.csv"
    noisy.to_csv(out_path)
    print(f"k={k}: {noisy.shape} -> {out_path}")

for omic in ["MET", "CNV", "OS"]:
    shutil.copyfile(f"{DATA_DIR}/{CANCER}_{omic}.csv",
                    f"{OUT_DIR}/{CANCER}_{omic}.csv")
    print(f"Copied {omic}")

print("Done.")
