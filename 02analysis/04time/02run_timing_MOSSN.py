"""Time MOSSN for one synthetic timing grid point."""

import sys

import pandas as pd

sys.path.insert(0, "/proj/c.zihao/work1/function/ablation")
from timing_common import run_timing
import MOSSN_func as MOSSN

def run_mossn(expr_file, links_file):
    links = pd.read_csv(links_file, index_col=0)
    expression = pd.read_csv(expr_file, index_col=0)
    graph, weights, expression = MOSSN.MOSSN_prepare_data(
        links=links, expression_data=expression
    )
    for sample_id in expression.columns:
        MOSSN.MOSSN_single_sample(
            sample_id=sample_id,
            G=graph,
            uniform_original_weights=weights,
            expression_data=expression,
        )

if __name__ == "__main__":
    run_timing("MOSSN", run_mossn)
