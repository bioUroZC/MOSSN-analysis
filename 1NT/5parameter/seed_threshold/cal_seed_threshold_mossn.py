from pathlib import Path

import pandas as pd

import sys
sys.path.append("/proj/c.zihao/work1/function/")
import MOSSN_noPrior


BASE_DIR = Path("/proj/c.zihao/work1/1NT/5parameter/seed_threshold")
DATA_FILE = Path("/proj/c.zihao/work1/1NT/5parameter/data/LUAD_paired_expr.csv")
LINK_FILE = Path("/proj/c.zihao/work1/1NT/1data/string/links.csv")

SEED_THRESHOLDS = [0.70, 0.80, 0.90, 0.95]
RESTART_ALPHA = 0.30
GAMMA = 2.0


def sanitize_threshold(threshold: float) -> str:
    return f"q{int(round(threshold * 100)):02d}"


def build_matrix(expression_data, links, alpha, seed_threshold):
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
            gamma=GAMMA,
            rwr_alpha=alpha,
            seed_quantile=seed_threshold
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

    for seed_threshold in SEED_THRESHOLDS:
        threshold_dir = BASE_DIR / sanitize_threshold(seed_threshold)
        threshold_dir.mkdir(parents=True, exist_ok=True)

        print(
            f"[seed-threshold] MOSSN_noPrior "
            f"alpha={RESTART_ALPHA} quantile={seed_threshold}"
        )
        matrix_df = build_matrix(
            expression_data=expression_data.copy(),
            links=links.copy(),
            alpha=RESTART_ALPHA,
            seed_threshold=seed_threshold
        )
        matrix_df.to_csv(threshold_dir / "merged_matrix.csv")


if __name__ == "__main__":
    main()
