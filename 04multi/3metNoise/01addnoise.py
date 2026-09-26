"""
Generate noisy MET data for LUAD at multiple noise levels.

Only the MET matrix is perturbed:
    noise_std = k * per-gene std(MET across samples)
Values are clipped to [0, 1], the valid beta-value range, so that the
perturbed data stay physically interpretable.

EXP, CNV, and OS are copied unchanged into 01data/ for convenience.
"""
import os, shutil
import numpy as np
import pandas as pd

CANCER      = "LUAD"
DATA_DIR    = "/proj/c.zihao/work1/04multi/1prepare/06files"
OUT_DIR     = "/proj/c.zihao/work1/04multi/3metNoise/01data"
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
SEED        = 42
CLIP        = (0.0, 1.0)

os.makedirs(OUT_DIR, exist_ok=True)

met = pd.read_csv(f"{DATA_DIR}/{CANCER}_MET.csv", index_col=0)
gene_stds = met.std(axis=1)
rng = np.random.default_rng(SEED)

for k in NOISE_LEVELS:
    if k == 0.0:
        noisy = met.copy()
    else:
        noise = rng.normal(0, 1, size=met.shape) * (gene_stds.values[:, None] * k)
        noisy = (met + noise).clip(lower=CLIP[0], upper=CLIP[1])
    out_path = f"{OUT_DIR}/{CANCER}_MET_k{k}.csv"
    noisy.to_csv(out_path)
    print(f"k={k}: {noisy.shape} -> {out_path}")

for omic in ["EXP", "CNV", "OS"]:
    shutil.copyfile(f"{DATA_DIR}/{CANCER}_{omic}.csv",
                    f"{OUT_DIR}/{CANCER}_{omic}.csv")
    print(f"Copied {omic}")

print("Done.")
