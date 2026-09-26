"""
Generate noisy EXP, MET, and CNV data for LUAD at multiple noise levels.

All three layers are perturbed simultaneously:
    noise_std = k * per-gene std(layer across samples)
MET is clipped to [0, 1] (beta values) and CNV to [0, 4] (copy number), so the
perturbed data stay physically interpretable.

Each layer draws from its own RNG stream seeded with SEED, in the same order as
the single-layer pipelines (3expNoise, 3metNoise, 3cnvNoise). The noise applied
to a given layer at a given k is therefore identical across all four pipelines,
so the single-layer and all-layer experiments differ only in which layers are
perturbed.

OS is copied unchanged into 01data/ for convenience.
"""
import os, shutil
import numpy as np
import pandas as pd

CANCER      = "LUAD"
DATA_DIR    = "/proj/c.zihao/work1/04multi/1prepare/06files"
OUT_DIR     = "/proj/c.zihao/work1/04multi/4allNoise/01data"
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
SEED        = 42
CLIP        = {"EXP": None, "MET": (0.0, 1.0), "CNV": (0.0, 4.0)}

os.makedirs(OUT_DIR, exist_ok=True)

for omic in ["EXP", "MET", "CNV"]:
    data = pd.read_csv(f"{DATA_DIR}/{CANCER}_{omic}.csv", index_col=0)
    gene_stds = data.std(axis=1)
    rng = np.random.default_rng(SEED)
    bounds = CLIP[omic]

    for k in NOISE_LEVELS:
        if k == 0.0:
            noisy = data.copy()
        else:
            noise = rng.normal(0, 1, size=data.shape) * (gene_stds.values[:, None] * k)
            noisy = data + noise
            if bounds is not None:
                noisy = noisy.clip(lower=bounds[0], upper=bounds[1])
        out_path = f"{OUT_DIR}/{CANCER}_{omic}_k{k}.csv"
        noisy.to_csv(out_path)
        print(f"{omic} k={k}: {noisy.shape} -> {out_path}")

shutil.copyfile(f"{DATA_DIR}/{CANCER}_OS.csv", f"{OUT_DIR}/{CANCER}_OS.csv")
print("Copied OS")
print("Done.")
