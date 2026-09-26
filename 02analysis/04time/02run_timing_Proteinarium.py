"""Time Proteinarium for one synthetic timing grid point."""

import os
import shutil
import sys
import tempfile

sys.path.insert(0, "/proj/c.zihao/work1/function/benchmark")
from timing_common import run_timing
import Proteinarium

TMPFS = "/dev/shm" if os.path.isdir("/dev/shm") else "/tmp"

def run_proteinarium(expr_file, links_file):
    output_dir = tempfile.mkdtemp(dir=TMPFS, prefix="t_Proteinarium_")
    try:
        Proteinarium.ProteinariumCal(
            exprSetFile=expr_file, link_file=links_file, save_path=output_dir
        )
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)

if __name__ == "__main__":
    run_timing("Proteinarium", run_proteinarium)
