"""Time SSN only, using the TCGA-LUAD paired-normal reference.

Usage:
  python 02run_timing_paired_ssn.py N_SAMPLES N_EDGES

The synthetic expression and edge matrices are reused from ``data/``.  Results
are written separately so the existing timing measurements for methods that do
not use a normal reference are preserved.
"""

import csv
import importlib
import os
import shutil
import sys
import tempfile
import time
import tracemalloc

if len(sys.argv) != 3:
    print("Usage: python 02run_timing_paired_ssn.py N_SAMPLES N_EDGES")
    sys.exit(1)

N_SAMPLES = int(sys.argv[1])
N_EDGES = int(sys.argv[2])

sys.path.insert(0, "/proj/c.zihao/work1/function/benchmark")

BASE_DIR = "/proj/c.zihao/work1/02analysis/04time"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
NORMAL_REF = "/proj/c.zihao/work1/01data/03TCGAnormal/paired_normals/TCGA-LUAD_paired_normal.csv"
TMPFS = "/dev/shm" if os.path.isdir("/dev/shm") else "/tmp"
N_REPS = 3

EXPR_FILE = os.path.join(DATA_DIR, f"expr_S{N_SAMPLES:03d}.csv")
LINKS_FILE = os.path.join(DATA_DIR, f"links_E{N_EDGES:06d}.csv")
OUT_CSV = os.path.join(
    RESULTS_DIR, f"timing_pairedSSN_S{N_SAMPLES:03d}_E{N_EDGES:06d}.csv"
)
FIELDNAMES = [
    "method", "n_samples", "n_edges", "rep",
    "wall_time_s", "time_per_sample_s", "peak_memory_mb",
]

def _tmpdir():
    return tempfile.mkdtemp(dir=TMPFS, prefix="t_pairedSSN_")

def _measure(function):
    tracemalloc.start()
    start = time.perf_counter()
    function()
    wall_time = time.perf_counter() - start
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return wall_time, peak_bytes / 1024 / 1024

def _write_row(row):
    exists = os.path.isfile(OUT_CSV)
    with open(OUT_CSV, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerow(row)

def run_ssn():
    import SSN

    output_dir = _tmpdir()
    try:
        return _measure(lambda: SSN.SSNcal(
            exprSetFile=EXPR_FILE,
            NormalFile=NORMAL_REF,
            link_file=LINKS_FILE,
            organ=None,
            save_path=output_dir,
        ))
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for path in (EXPR_FILE, LINKS_FILE, NORMAL_REF):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Required input not found: {path}")

    print(f"N_SAMPLES={N_SAMPLES} N_EDGES={N_EDGES} -> {OUT_CSV}", flush=True)
    if os.path.isfile(OUT_CSV):
        os.remove(OUT_CSV)

    importlib.import_module("SSN")
    for rep in range(N_REPS):
        print(f"  [SSN, paired normal] rep={rep} ...", flush=True)
        try:
            wall_time, peak_mb = run_ssn()
            _write_row({
                "method": "SSN",
                "n_samples": N_SAMPLES,
                "n_edges": N_EDGES,
                "rep": rep,
                "wall_time_s": round(wall_time, 4),
                "time_per_sample_s": round(wall_time / N_SAMPLES, 6),
                "peak_memory_mb": round(peak_mb, 2),
            })
            print(
                f"    wall={wall_time:.2f}s /sample={wall_time / N_SAMPLES:.4f}s "
                f"peak={peak_mb:.1f}MB",
                flush=True,
            )
        except Exception as error:
            print(f"    ERROR: {error}", flush=True)
            _write_row({
                "method": "SSN",
                "n_samples": N_SAMPLES,
                "n_edges": N_EDGES,
                "rep": rep,
                "wall_time_s": "ERROR",
                "time_per_sample_s": "ERROR",
                "peak_memory_mb": "ERROR",
            })

    print(f"Done -> {OUT_CSV}")

if __name__ == "__main__":
    main()
