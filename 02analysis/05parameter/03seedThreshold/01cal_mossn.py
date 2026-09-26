from pathlib import Path
import sys

import pandas as pd

sys.path.append("/proj/c.zihao/work1/function/ablation")
import MOSSN_func as MOSSN

BASE_DIR = Path("/proj/c.zihao/work1/02analysis/05parameter/03seedThreshold")
DATA_FILE = Path("/proj/c.zihao/work1/02analysis/05parameter/00data/LUAD_paired_expr.csv")
LINK_FILE = Path("/proj/c.zihao/work1/01data/06string/links.csv")
SEED_QUANTILES = [0.70, 0.80, 0.90, 0.95]
LAM = 2.0
RESTART_ALPHA = 0.30

def label_quantile(quantile: float) -> str:
    return f"q{int(round(quantile * 100)):02d}"

def build_matrix(expression_data: pd.DataFrame, links: pd.DataFrame, seed_quantile: float) -> pd.DataFrame:
    graph, uniform_weights, expression_data = MOSSN.MOSSN_prepare_data(links, expression_data)
    sample_columns = {}
    edge_index = None
    for sample_id in expression_data.columns:
        edge_df = MOSSN.MOSSN_single_sample(
            sample_id=sample_id, G=graph, uniform_original_weights=uniform_weights,
            expression_data=expression_data, lam=LAM, rwr_alpha=RESTART_ALPHA,
            seed_quantile=seed_quantile,
        )
        if edge_index is None:
            edge_index = edge_df["Node1"] + "_" + edge_df["Node2"]
        sample_columns[sample_id] = edge_df["FinalWeight"].to_numpy()
    return pd.DataFrame(sample_columns, index=edge_index)

def main() -> None:
    expression_data = pd.read_csv(DATA_FILE, index_col=0)
    links = pd.read_csv(LINK_FILE, index_col=0)
    for seed_quantile in SEED_QUANTILES:
        output_dir = BASE_DIR / label_quantile(seed_quantile)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"[seed threshold] MOSSN seed_quantile={seed_quantile}")
        build_matrix(expression_data.copy(), links.copy(), seed_quantile).to_csv(output_dir / "merged_matrix.csv")

if __name__ == "__main__":
    main()
