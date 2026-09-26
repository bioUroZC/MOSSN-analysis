"""Shared timing and CSV-writing helpers for the Python benchmark methods."""

import csv
import os
import sys
import time
import tracemalloc

BASE_DIR = "/proj/c.zihao/work1/02analysis/04time"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
N_REPS = 3
FIELDNAMES = [
    "method", "n_samples", "n_edges", "rep",
    "wall_time_s", "time_per_sample_s", "peak_memory_mb",
]

def run_timing(method, runner):
    """Time one method for a requested (n_samples, n_edges) grid point.

    ``runner`` receives ``(expr_file, links_file)`` and performs one complete
    method calculation.  It is called three times and the resulting records are
    written to a method-specific CSV so independent jobs cannot overwrite one
    another's output.
    """
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(sys.argv[0])} N_SAMPLES N_EDGES")
        sys.exit(1)

    n_samples = int(sys.argv[1])
    n_edges = int(sys.argv[2])
    expr_file = os.path.join(DATA_DIR, f"expr_S{n_samples:03d}.csv")
    links_file = os.path.join(DATA_DIR, f"links_E{n_edges:06d}.csv")
    out_csv = os.path.join(
        RESULTS_DIR, f"timing_{method}_S{n_samples:03d}_E{n_edges:06d}.csv"
    )

    for path in (expr_file, links_file):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Required input not found: {path}")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    if os.path.isfile(out_csv):
        os.remove(out_csv)

    def write_row(row):
        exists = os.path.isfile(out_csv)
        with open(out_csv, "a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
            if not exists:
                writer.writeheader()
            writer.writerow(row)

    print(f"{method}: N_SAMPLES={n_samples} N_EDGES={n_edges} -> {out_csv}", flush=True)
    for rep in range(N_REPS):
        print(f"  [{method}] rep={rep} ...", flush=True)
        try:
            tracemalloc.start()
            start = time.perf_counter()
            runner(expr_file, links_file)
            wall_time = time.perf_counter() - start
            _, peak_bytes = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            peak_mb = peak_bytes / 1024 / 1024
            write_row({
                "method": method,
                "n_samples": n_samples,
                "n_edges": n_edges,
                "rep": rep,
                "wall_time_s": round(wall_time, 4),
                "time_per_sample_s": round(wall_time / n_samples, 6),
                "peak_memory_mb": round(peak_mb, 2),
            })
            print(
                f"    wall={wall_time:.2f}s /sample={wall_time / n_samples:.4f}s "
                f"peak={peak_mb:.1f}MB",
                flush=True,
            )
        except Exception as error:
            if tracemalloc.is_tracing():
                tracemalloc.stop()
            print(f"    ERROR: {error}", flush=True)
            write_row({
                "method": method,
                "n_samples": n_samples,
                "n_edges": n_edges,
                "rep": rep,
                "wall_time_s": "ERROR",
                "time_per_sample_s": "ERROR",
                "peak_memory_mb": "ERROR",
            })

    print(f"Done -> {out_csv}")
