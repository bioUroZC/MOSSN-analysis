from pathlib import Path

import pandas as pd

import sys
sys.path.append("/proj/c.zihao/work1/function/")
import MOSSN_noPrior


BASE_DIR = Path("/proj/c.zihao/work1/1NT/5parameter/gamma")
DATA_FILE = Path("/proj/c.zihao/work1/1NT/5parameter/data/LUAD_paired_expr.csv")
LINK_FILE = Path("/proj/c.zihao/work1/1NT/1data/string/links.csv")

GAMMAS = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
RESTART_ALPHA = 0.30
SEED_THRESHOLD = 0.90


def sanitize_gamma(gamma: float) -> str:
    return f"g{int(round(gamma * 10)):03d}"


def build_matrix(expression_data, links, gamma):
    G, base_weights, expression_data = MOSSN_noPrior.prepare_data_MOSSN_noPrior(
        links=links,
        expression_data=expression_data
    )

    sample_columns = {}
    edge_index = None

    for sample_id in expression_data.columns:
        edge_df = MOSSN_noPrior.MOSSN_noPrior_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=base_weights,
            expression_data=expression_data,
            gamma=gamma,
            rwr_alpha=RESTART_ALPHA,
            seed_quantile=SEED_THRESHOLD
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

    for gamma in GAMMAS:
        gamma_dir = BASE_DIR / sanitize_gamma(gamma)
        gamma_dir.mkdir(parents=True, exist_ok=True)

        print(f"[gamma] MOSSN_noPrior gamma={gamma}")
        matrix_df = build_matrix(
            expression_data=expression_data.copy(),
            links=links.copy(),
            gamma=gamma
        )
        matrix_df.to_csv(gamma_dir / "merged_matrix.csv")


if __name__ == "__main__":
    main()
