import os
from pathlib import Path

import pandas as pd

import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import MOSSN_uniform


BASE_DIR = Path(PROJ_ROOT + "/1NT/5parameter/restart")
DATA_FILE = Path(PROJ_ROOT + "/1NT/5parameter/data/LUAD_paired_expr.csv")
LINK_FILE = Path(PROJ_ROOT + "/1NT/1data/string/links.csv")

ALPHAS = [round(x / 10, 1) for x in range(1, 10)]
GAMMA = 2.0
SEED_QUANTILE = 0.9


def sanitize_alpha(alpha: float) -> str:
    return f"alpha{int(round(alpha * 10)):02d}"


def build_matrix(expression_data, links, alpha):
    G, real_original_weights, expression_data = MOSSN_uniform.prepare_data_MOSSN_uniform(
        links=links,
        expression_data=expression_data
    )

    sample_columns = {}
    edge_index = None

    for sample_id in expression_data.columns:
        edge_df = MOSSN_uniform.MOSSN_uniform_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            gamma=GAMMA,
            rwr_alpha=alpha,
            seed_quantile=SEED_QUANTILE
        )
        edge_ids = edge_df["Node1"] + "_" + edge_df["Node2"]
        if edge_index is None:
            edge_index = edge_ids
        sample_columns[sample_id] = edge_df["FinalWeight"].to_numpy()

    return pd.DataFrame(sample_columns, index=edge_index)


def main():
    BASE_DIR.mkdir(parents=True, exist_ok=True)

    expression_data = pd.read_csv(DATA_FILE, index_col=0)
    links = pd.read_csv(LINK_FILE, index_col=0)

    for alpha in ALPHAS:
        alpha_dir = BASE_DIR / sanitize_alpha(alpha)
        alpha_dir.mkdir(parents=True, exist_ok=True)

        print(f"[restart] MOSSN_uniform alpha={alpha}")
        matrix_df = build_matrix(expression_data.copy(), links.copy(), alpha)
        matrix_df.to_csv(alpha_dir / "merged_matrix.csv")


if __name__ == "__main__":
    main()
