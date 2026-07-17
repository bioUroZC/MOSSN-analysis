#!/usr/bin/env python3
"""Single-point timing job: run all methods for one (n_samples, n_edges) pair.

Usage:
  python run_timing.py N_SAMPLES N_EDGES

Examples:
  python run_timing.py 10 10000    # sample dim, N=10
  python run_timing.py 50 10000    # sample dim, N=50
  python run_timing.py 10 50000    # network dim, E=50k

Output: results/timing_S{N:03d}_E{K:06d}.csv
"""

import os, sys, csv, time, shutil, tempfile, tracemalloc, importlib
import pandas as pd

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")


if len(sys.argv) != 3:
    print("Usage: python run_timing.py N_SAMPLES N_EDGES")
    sys.exit(1)

N_SAMPLES = int(sys.argv[1])
N_EDGES   = int(sys.argv[2])

sys.path.insert(0, PROJ_ROOT + "/function")

BASE_DIR    = PROJ_ROOT + "/1NT/6time"
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
NORMAL_REF  = PROJ_ROOT + "/0ref/Test8/combined_expr_df.csv"
ORGAN       = "Lung"
TMPFS       = "/dev/shm" if os.path.isdir("/dev/shm") else "/tmp"
N_REPS      = 3

os.makedirs(RESULTS_DIR, exist_ok=True)

EXPR_FILE    = os.path.join(DATA_DIR, f"expr_S{N_SAMPLES:03d}.csv")
LINKS_FILE   = os.path.join(DATA_DIR, f"links_E{N_EDGES:06d}.csv")
OUT_CSV      = os.path.join(RESULTS_DIR, f"timing_S{N_SAMPLES:03d}_E{N_EDGES:06d}.csv")

FIELDNAMES = ["method", "n_samples", "n_edges", "rep",
              "wall_time_s", "time_per_sample_s", "peak_memory_mb"]

# ── helpers ────────────────────────────────────────────────────────────────────

def _tmpdir(tag):
    return tempfile.mkdtemp(dir=TMPFS, prefix=f"t_{tag}_")

def _measure(fn):
    tracemalloc.start()
    t0 = time.perf_counter()
    fn()
    wall = time.perf_counter() - t0
    _, peak_b = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return wall, peak_b / 1024 / 1024

def _write_row(row):
    exists = os.path.isfile(OUT_CSV)
    with open(OUT_CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if not exists:
            w.writeheader()
        w.writerow(row)

# ── runners ────────────────────────────────────────────────────────────────────

def run_SWEET():
    import SWEET
    d = _tmpdir("SWEET")
    try:
        return _measure(lambda: SWEET.SWEETcal(
            exprSetFile=EXPR_FILE, link_file=LINKS_FILE, save_path=d))
    finally:
        shutil.rmtree(d, ignore_errors=True)

def run_SSN():
    import SSN
    d = _tmpdir("SSN")
    try:
        return _measure(lambda: SSN.SSNcal(
            exprSetFile=EXPR_FILE, NormalFile=NORMAL_REF,
            link_file=LINKS_FILE, organ=ORGAN, save_path=d))
    finally:
        shutil.rmtree(d, ignore_errors=True)

def run_Patkar():
    import Patkar
    d = _tmpdir("Patkar")
    try:
        return _measure(lambda: Patkar.PatkarCal(
            exprSetFile=EXPR_FILE, link_file=LINKS_FILE, save_path=d))
    finally:
        shutil.rmtree(d, ignore_errors=True)

def run_PPIXpress():
    import PPIXpress
    d = _tmpdir("PPIXpress")
    try:
        return _measure(lambda: PPIXpress.PPIXpressCal(
            exprSetFile=EXPR_FILE, link_file=LINKS_FILE, save_path=d))
    finally:
        shutil.rmtree(d, ignore_errors=True)

def run_Proteinarium():
    import Proteinarium
    d = _tmpdir("Proteinarium")
    try:
        return _measure(lambda: Proteinarium.ProteinariumCal(
            exprSetFile=EXPR_FILE, link_file=LINKS_FILE, save_path=d))
    finally:
        shutil.rmtree(d, ignore_errors=True)

def run_MOSSN_uniform():
    import MOSSN_uniform  # resolved at runtime via sys.path
    def _run():
        links = pd.read_csv(LINKS_FILE,   index_col=0)
        expr  = pd.read_csv(EXPR_FILE,    index_col=0)
        G, weights, e = MOSSN_uniform.prepare_data_MOSSN_uniform(
            links=links, expression_data=expr)
        for sid in e.columns:
            MOSSN_uniform.MOSSN_uniform_single_sample(
                sample_id=sid, G=G, real_original_weights=weights,
                expression_data=e)
    return _measure(_run)

METHODS = {
    "SWEET":        run_SWEET,
    "SSN":          run_SSN,
    "Patkar":       run_Patkar,
    "PPIXpress":    run_PPIXpress,
    "Proteinarium": run_Proteinarium,
    "MOSSN_uniform": run_MOSSN_uniform,
}

MODULES = {
    "SWEET": "SWEET",
    "SSN": "SSN",
    "Patkar": "Patkar",
    "PPIXpress": "PPIXpress",
    "Proteinarium": "Proteinarium",
    "MOSSN_uniform": "MOSSN_uniform",
}

# ── main ───────────────────────────────────────────────────────────────────────

print(f"N_SAMPLES={N_SAMPLES}  N_EDGES={N_EDGES}  → {OUT_CSV}", flush=True)

# Replace any previous partial or duplicate output when rerunning this grid point.
if os.path.isfile(OUT_CSV):
    os.remove(OUT_CSV)

# Import once before timing so the measured region focuses on runtime rather than
# one-time module import overhead.
for method_name, runner in METHODS.items():
    print(f"  warming import for [{method_name}] ...", flush=True)
    try:
        importlib.import_module(MODULES[method_name])
    except Exception as e:
        print(f"    warm-up skipped due to: {e}", flush=True)

for method_name, runner in METHODS.items():
    for rep in range(N_REPS):
        print(f"  [{method_name}] rep={rep} ...", flush=True)
        try:
            wall, peak_mb = runner()
            row = dict(method=method_name, n_samples=N_SAMPLES, n_edges=N_EDGES,
                       rep=rep, wall_time_s=round(wall, 4),
                       time_per_sample_s=round(wall / N_SAMPLES, 6),
                       peak_memory_mb=round(peak_mb, 2))
            _write_row(row)
            print(f"    wall={wall:.2f}s  /sample={wall/N_SAMPLES:.4f}s  peak={peak_mb:.1f}MB")
        except Exception as e:
            print(f"    ERROR: {e}")
            _write_row(dict(method=method_name, n_samples=N_SAMPLES, n_edges=N_EDGES,
                            rep=rep, wall_time_s="ERROR",
                            time_per_sample_s="ERROR", peak_memory_mb="ERROR"))

print(f"Done → {OUT_CSV}")
