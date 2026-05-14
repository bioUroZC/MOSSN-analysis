#!/usr/bin/env python3
"""Merge per-grid timing CSVs into plot-ready result tables.

Outputs:
  results/timing_results_py.csv
  results/timing_results_LIONESS.csv
"""

import glob
import os
import pandas as pd


BASE_DIR = "/proj/c.zihao/work1/1NT/6time"
RESULTS_DIR = os.path.join(BASE_DIR, "results")


def _concat(pattern):
    files = sorted(glob.glob(os.path.join(RESULTS_DIR, pattern)))
    if not files:
        raise FileNotFoundError(f"No result files matched: {pattern}")
    frames = [pd.read_csv(path) for path in files]
    return pd.concat(frames, ignore_index=True)


def main():
    py = _concat("timing_S*_E*.csv")
    li = _concat("timing_LIONESS_S*_E*.csv")

    py.to_csv(os.path.join(RESULTS_DIR, "timing_results_py.csv"), index=False)
    li.to_csv(os.path.join(RESULTS_DIR, "timing_results_LIONESS.csv"), index=False)

    print("Saved:")
    print(f"  {os.path.join(RESULTS_DIR, 'timing_results_py.csv')}")
    print(f"  {os.path.join(RESULTS_DIR, 'timing_results_LIONESS.csv')}")


if __name__ == "__main__":
    main()
