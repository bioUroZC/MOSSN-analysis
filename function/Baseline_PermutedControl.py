import hashlib

import numpy as np


def seed_from_sample_id(sample_id, base_seed=1):
    """Deterministic but sample-varying seed. Each sample must get an
    independent shuffle - reusing one fixed seed across every sample would
    apply the exact same permutation index to every sample, which merely
    relabels edges (each edge's cross-sample value vector is preserved
    under a new name) instead of actually scrambling the sample-specific
    edge <-> weight correspondence. Uses hashlib instead of Python's
    built-in hash() so the seed is reproducible across runs/processes
    (str hashing is randomized per-process unless PYTHONHASHSEED is set).
    """
    digest = hashlib.md5(f"{base_seed}_{sample_id}".encode()).hexdigest()
    return int(digest, 16) % (2**32 - 1)


def permute_single_sample_edges(edge_df, seed=1, value_col="FinalWeight"):
    """Post-hoc permuted control: given one sample's already-computed full
    MOSSN edge table, randomly permute which edge receives which weight
    value. Topology and the sample's score distribution are unchanged;
    only the edge <-> score assignment is scrambled. No RWR/correction
    needs to be rerun.

    `seed` must vary per sample (see seed_from_sample_id) - a shared seed
    across samples would make this a no-op control (see docstring above).
    """
    rng = np.random.default_rng(seed)
    permuted = edge_df.copy()
    permuted[value_col] = rng.permutation(permuted[value_col].to_numpy())
    return permuted
