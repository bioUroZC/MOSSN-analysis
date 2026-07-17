import numpy as np
import pandas as pd


def EdgeNoRWR_single_sample(sample_id, G, uniform_original_weights, expression_data, gamma=2.0):
    for (u, v), weight in uniform_original_weights.items():
        if G.has_edge(u, v):
            G[u][v]["weight"] = weight
            G[u][v]["adjusted_weight"] = weight
            G[u][v]["correction_score"] = 1.0

    sample_expression = expression_data[sample_id]
    sample_median = sample_expression.median()
    sample_iqr = sample_expression.quantile(0.75) - sample_expression.quantile(0.25)
    if sample_iqr == 0:
        sample_iqr = 1e-6

    expr_norm = ((sample_expression - sample_median) / sample_iqr).copy()
    alpha_mod = 1.0

    for u, v, data in G.edges(data=True):
        zu = expr_norm.get(u, 0)
        zv = expr_norm.get(v, 0)
        m_uv = zu + zv
        sigmoid = 1 / (1 + np.exp(-alpha_mod * m_uv))
        correction_score = 1 + gamma * (sigmoid - 0.5)
        correction_score = max(correction_score, 0)
        base_weight = uniform_original_weights.get((u, v), uniform_original_weights.get((v, u), 0))
        data["adjusted_weight"] = base_weight * correction_score
        data["correction_score"] = correction_score

    # Ablation: no RWR -> final weight is the corrected edge weight itself.
    edges_data = []
    for u, v, data in G.edges(data=True):
        base_weight = uniform_original_weights.get((u, v), uniform_original_weights.get((v, u), 0))
        edges_data.append(
            {
                "Sample": sample_id,
                "Node1": u,
                "Node2": v,
                "BaseWeight": base_weight,
                "FinalWeight": data["adjusted_weight"],
            }
        )

    edge_weights_df = pd.DataFrame(edges_data)
    edge_weights_df["FinalWeight"] = edge_weights_df["FinalWeight"].round(5)
    return edge_weights_df
