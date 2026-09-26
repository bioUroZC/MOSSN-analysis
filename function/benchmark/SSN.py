import numpy as np
import pandas as pd
from scipy.stats import pearsonr, norm
from collections import defaultdict
from tqdm import tqdm

def SSNcal(exprSetFile, NormalFile, link_file, organ=None, save_path=None):
    """Calculate SSN edge scores from a normal reference expression matrix.

    Parameters
    ----------
    NormalFile
        Either (i) a gene-by-sample normal-expression matrix, or (ii) the
        legacy GTEx combined matrix containing an ``organ`` column.
    organ
        Required only for the legacy GTEx combined matrix.  Set to ``None``
        when ``NormalFile`` is already the relevant gene-by-sample reference
        matrix (e.g. TCGA paired adjacent-normal samples).
    """
    if save_path is None:
        raise ValueError("save_path must be provided.")

    Normal = pd.read_csv(NormalFile, index_col=0)
    if organ is None:
        
        if "organ" in Normal.columns:
            raise ValueError(
                "organ=None requires a gene-by-sample normal matrix without "
                "an 'organ' column. Supply organ=<tissue> for the legacy GTEx file."
            )
        ref_df = Normal
        reference_label = "provided normal matrix"
    else:
        
        if "organ" not in Normal.columns:
            raise ValueError(
                "The supplied normal matrix has no 'organ' column; set organ=None "
                "to use it directly."
            )
        organ_df = Normal[Normal["organ"] == organ]
        organ_df = organ_df.drop(columns=["organ"]).T
        ref_df = organ_df
        reference_label = f"GTEx {organ}"

    n_ref_samples = ref_df.shape[1]
    if n_ref_samples < 2:
        raise ValueError(
            f"SSN requires at least 2 normal reference samples for {reference_label}, "
            f"but got {n_ref_samples} from {NormalFile}. "
            "Increase the normal subset size before running SSN."
        )

    initial_gene_count = ref_df.shape[0]
    ref_df = ref_df.loc[ref_df.std(axis=1) > 1e-5]
    filtered_gene_count = ref_df.shape[0]
    print(f"Filtered constant genes: {initial_gene_count - filtered_gene_count} removed, {filtered_gene_count} remaining")

    if filtered_gene_count == 0:
        raise ValueError(
            f"No reference genes remained after variance filtering for organ '{organ}'. "
            f"Check whether {NormalFile} contains enough non-constant normal samples."
        )

    test_df = pd.read_csv(exprSetFile, index_col=0)
    test_df = test_df.loc[test_df.std(axis=1) > 0]

    common_genes = ref_df.index.intersection(test_df.index)

    links = pd.read_csv(link_file, index_col=0)
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = common_genes.intersection(genes_in_links)

    links = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
        ].drop_duplicates(subset=['protein1', 'protein2'])

    print(f"Filtered STRING links: {links.shape}")

    used_genes = pd.unique(links[['protein1', 'protein2']].values.ravel())
    ref_aligned = ref_df.loc[used_genes]
    test_aligned = test_df.loc[used_genes]
    print(f"Number of genes used in expression: {len(used_genes)}")

    gene_pairs = [
        (g1, g2)
        for g1, g2 in zip(links["protein1"], links["protein2"])
        if g1 in ref_aligned.index and g2 in ref_aligned.index
    ]
    print(f"Total STRING gene pairs with expression data: {len(gene_pairs)}")

    def compute_reference_pcc(reference_df, description):
        """Calculate one baseline PCC value per STRING edge."""
        pcc_ref = {}
        for g1, g2 in tqdm(gene_pairs, desc=description):
            try:
                r, _ = pearsonr(reference_df.loc[g1], reference_df.loc[g2])
            except Exception:
                r = 0
            pcc_ref[(g1, g2)] = r
        return pcc_ref

    reference_sample_ids = set(ref_aligned.columns)
    test_sample_ids = list(test_aligned.columns)
    full_reference_pcc = None
    if any(sample not in reference_sample_ids for sample in test_sample_ids):
        print("Computing full-reference PCCs for non-reference test samples...")
        full_reference_pcc = compute_reference_pcc(ref_aligned, "Full reference PCC")

    def compute_sample_stat(sample_name):
        if sample_name in reference_sample_ids:
            sample_reference = ref_aligned.drop(columns=sample_name)
            n_reference = sample_reference.shape[1]
            if n_reference < 2:
                raise ValueError(
                    f"Leave-one-out SSN requires at least 3 normal samples, but "
                    f"'{sample_name}' leaves only {n_reference} in {NormalFile}."
                )
            pcc_reference = compute_reference_pcc(
                sample_reference,
                f"LOO reference PCC: {sample_name}"
            )
        else:
            sample_reference = ref_aligned
            n_reference = sample_reference.shape[1]
            pcc_reference = full_reference_pcc

        sample_df = test_aligned[[sample_name]]
        combined_df = pd.concat([sample_reference, sample_df], axis=1)

        z_scores, p_values, delta_pccs = {}, {}, {}

        for edge, pcc_n in pcc_reference.items():
            g1, g2 = edge
            try:
                r, _ = pearsonr(combined_df.loc[g1], combined_df.loc[g2])
            except Exception:
                r = 0
            delta = r - pcc_n
            
            sigma = (1 - pcc_n ** 2) / np.sqrt(n_reference - 1)
            if sigma == 0:
                continue
            z = delta / sigma
            p = 2 * (1 - norm.cdf(abs(z)))  

            z_scores[edge] = round(z, 3)
            delta_pccs[edge] = round(delta, 3)
            p_values[edge] = p

        return sample_name, z_scores, delta_pccs, p_values

    print("Computing stats for test samples (leave-one-out for reference normals)...")
    results = []
    for sample in tqdm(test_aligned.columns, desc="Samples"):
        results.append(compute_sample_stat(sample))

    z_scores_dict = defaultdict(dict)
    delta_pcc_dict = defaultdict(dict)
    p_value_dict = defaultdict(dict)

    for sample_name, z_scores, delta_pccs, p_values in results:
        for edge, z in z_scores.items():
            z_scores_dict[edge][sample_name] = z
        for edge, d in delta_pccs.items():
            delta_pcc_dict[edge][sample_name] = d
        for edge, p in p_values.items():
            p_value_dict[edge][sample_name] = p

    z_df = pd.DataFrame.from_dict(z_scores_dict, orient='index')
    z_df.index = pd.MultiIndex.from_tuples(z_df.index, names=["Gene1", "Gene2"])

    delta_df = pd.DataFrame.from_dict(delta_pcc_dict, orient='index')
    delta_df.index = pd.MultiIndex.from_tuples(delta_df.index, names=["Gene1", "Gene2"])

    pval_df = pd.DataFrame.from_dict(p_value_dict, orient='index')
    pval_df.index = pd.MultiIndex.from_tuples(pval_df.index, names=["Gene1", "Gene2"])

    z_df.to_csv(f"{save_path}/Zscore.csv")
    delta_df.to_csv(f"{save_path}/delta.csv")
    pval_df.to_csv(f"{save_path}/pvalue.csv")
