rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)

# Merges the two topology-ablation methods only: RandomBackbone,
# DegreePreservedBackbone. These rewire the PPI topology, so their edges
# mostly do not coincide with the real STRING backbone (`links_ref` in
# 3merge.R). Merging them against `links_ref` like the other edge methods
# would drop almost all of their edges as non-matches, so each is indexed
# by the union of interactions it actually generated instead.
#
# Run this standalone (the other methods in 3merge.R have already been
# merged and should not be re-run).

BASE_DIR <- paste0(PROJ_ROOT, "/1NT/2string/ablation")

available_datasets <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

topology_methods <- c("RandomBackbone", "DegreePreservedBackbone")

for (method in topology_methods) {
    all_df <- list()

    for (dataset in available_datasets) {
        folder <- file.path(BASE_DIR, method, dataset)
        files <- list.files(folder, pattern = "_edges\\.csv$", full.names = TRUE)
        if (length(files) == 0) next

        one_dataset <- lapply(files, readr::read_csv, show_col_types = FALSE) %>%
            dplyr::bind_rows()

        one_dataset <- one_dataset %>%
            dplyr::transmute(
                protein1 = Node1,
                protein2 = Node2,
                Sample = Sample,
                Weight = FinalWeight
            )

        all_df[[dataset]] <- one_dataset
    }

    method_long <- dplyr::bind_rows(all_df)
    rm(all_df)
    if (nrow(method_long) > 0) {
        sample_ids <- sort(unique(method_long$Sample))

        method_long <- method_long %>%
            dplyr::transmute(
                protein1_raw = protein1,
                protein2_raw = protein2,
                protein1 = pmin(protein1_raw, protein2_raw),
                protein2 = pmax(protein1_raw, protein2_raw),
                Sample = Sample,
                Weight = Weight
            )

        interaction_ids <- method_long %>%
            dplyr::transmute(
                Interaction = paste(protein1, protein2, sep = "_")
            ) %>%
            dplyr::distinct() %>%
            dplyr::pull(Interaction) %>%
            sort()

        out_mat <- matrix(
            0,
            nrow = length(interaction_ids),
            ncol = length(sample_ids)
        )
        rownames(out_mat) <- interaction_ids
        colnames(out_mat) <- sample_ids

        current_interaction <- paste(method_long$protein1, method_long$protein2, sep = "_")
        row_id <- match(current_interaction, interaction_ids)
        col_id <- match(method_long$Sample, sample_ids)
        out_mat[cbind(row_id, col_id)] <- method_long$Weight
        rm(method_long)

        method_matrix <- as.data.frame(out_mat, check.names = FALSE)
        method_matrix$Interaction <- rownames(out_mat)
        method_matrix <- method_matrix[, c("Interaction", sample_ids)]
        rownames(method_matrix) <- NULL

        readr::write_csv(method_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
        cat(method, "->", nrow(method_matrix), "interactions x",
            ncol(method_matrix) - 1, "samples\n")
        rm(method_matrix)
    } else {
        rm(method_long)
        cat(method, ": no input files found\n")
    }
    gc()
}
