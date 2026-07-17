rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

# 1. Basic settings
BASE_DIR <- paste0(PROJ_ROOT, "/1NT/2string/ablation")
LINKS_FILE <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")

available_datasets <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

# Edge-level ablation methods: per-sample "*_edges.csv"
#   (Sample, Node1, Node2, BaseWeight, FinalWeight)
# NOTE: RandomBackbone and DegreePreservedBackbone are handled separately
# below (they generate their own edge topology, which mostly does not
# coincide with the real STRING backbone, so they can't be indexed by
# `links_ref` like the methods below).
edge_methods <- c(
    "EdgeNoRWR", "MOSSN_noSeed", "MOSSN_noCorr", "PermutedControl"
)

# Node-level ablation methods: already a per-tumor wide matrix
#   "{method}/{dataset}.csv" with first column "Interaction" (feature) + sample columns
node_methods <- c("RawExpr", "NodeRWR")

# 2. Read reference interactions (fixed row set shared by all edge methods)
links_ref <- readr::read_csv(LINKS_FILE, show_col_types = FALSE) %>%
    dplyr::select(protein1, protein2) %>%
    dplyr::transmute(
        protein1_raw = protein1,
        protein2_raw = protein2,
        protein1 = pmin(protein1_raw, protein2_raw),
        protein2 = pmax(protein1_raw, protein2_raw)
    ) %>%
    dplyr::select(protein1, protein2) %>%
    dplyr::distinct()

links_ref$Interaction <- paste(links_ref$protein1, links_ref$protein2, sep = "_")

# 3. Helper: convert an edge long table to interaction x sample matrix
make_matrix_from_long <- function(df_long) {
    sample_ids <- sort(unique(df_long$Sample))
    interaction_ids <- links_ref$Interaction

    df_long <- df_long %>%
        dplyr::transmute(
            protein1_raw = protein1,
            protein2_raw = protein2,
            protein1 = pmin(protein1_raw, protein2_raw),
            protein2 = pmax(protein1_raw, protein2_raw),
            Sample = Sample,
            Weight = Weight
        )

    out_mat <- matrix(
        0,
        nrow = length(interaction_ids),
        ncol = length(sample_ids)
    )

    rownames(out_mat) <- interaction_ids
    colnames(out_mat) <- sample_ids

    current_interaction <- paste(df_long$protein1, df_long$protein2, sep = "_")
    row_id <- match(current_interaction, interaction_ids)
    col_id <- match(df_long$Sample, sample_ids)
    keep <- !is.na(row_id) & !is.na(col_id)

    out_mat[cbind(row_id[keep], col_id[keep])] <- df_long$Weight[keep]

    out_df <- as.data.frame(out_mat, check.names = FALSE)
    out_df$Interaction <- rownames(out_mat)
    out_df <- out_df[, c("Interaction", sample_ids)]
    rownames(out_df) <- NULL
    out_df
}

# ===================================================
# Edge-level methods: merge all tumors -> one interaction x sample matrix
# ===================================================
for (method in edge_methods) {
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
        method_matrix <- make_matrix_from_long(method_long)
        rm(method_long)
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

# ===================================================
# Node-level methods: merge all tumors -> one feature x sample matrix
# ===================================================
for (method in node_methods) {
    all_df <- list()

    for (dataset in available_datasets) {
        file <- file.path(BASE_DIR, method, paste0(dataset, ".csv"))
        if (!file.exists(file)) next

        one_dataset <- readr::read_csv(file, show_col_types = FALSE)
        names(one_dataset)[1] <- "Interaction"
        all_df[[dataset]] <- one_dataset
    }

    if (length(all_df) > 0) {
        method_matrix <- purrr::reduce(
            all_df,
            dplyr::full_join,
            by = "Interaction"
        )
        rm(all_df)
        # Fill missing sample values (genes absent in some tumors) with 0
        method_matrix <- method_matrix %>%
            dplyr::mutate(dplyr::across(-Interaction, ~ tidyr::replace_na(.x, 0)))

        readr::write_csv(method_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
        cat(method, "->", nrow(method_matrix), "features x",
            ncol(method_matrix) - 1, "samples\n")
        rm(method_matrix)
    } else {
        rm(all_df)
        cat(method, ": no input files found\n")
    }
    gc()
}
