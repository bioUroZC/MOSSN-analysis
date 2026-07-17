def get_raw_expression_features(expression_data):
    """Raw gene expression as a sample-by-gene feature matrix - the
    simplest possible "no network" baseline. Returns the matrix unchanged
    (genes x samples); kept as a function for symmetry with the other
    baseline modules and as a single place to adjust if a transform is
    ever needed.
    """
    return expression_data.copy()
