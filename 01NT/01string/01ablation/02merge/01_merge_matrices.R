rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

WORKFLOW_DIR <- "/proj/c.zihao/work1/01NT/01string"
BASE_DIR <- file.path(WORKFLOW_DIR, "01ablation")
LINKS_FILE <- "/proj/c.zihao/work1/01data/06string/links.csv"
OUTPUT_DIR <- file.path(BASE_DIR, "02merge", "merged_matrices")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

available_datasets <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

edge_methods <- c(
    "noRWR", "noSeed", "noCorr", "PermutedControl"
)

node_methods <- c("RawExpr", "NodeRWR")

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
        method_output_dir <- file.path(OUTPUT_DIR, method)
        dir.create(method_output_dir, recursive = TRUE, showWarnings = FALSE)
        readr::write_csv(method_matrix, file.path(method_output_dir, "merged_matrix.csv"))
        cat(method, "->", nrow(method_matrix), "interactions x",
            ncol(method_matrix) - 1, "samples\n")
        rm(method_matrix)
    } else {
        rm(method_long)
        cat(method, ": no input files found\n")
    }
    gc()
}

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
        
        method_matrix <- method_matrix %>%
            dplyr::mutate(dplyr::across(-Interaction, ~ tidyr::replace_na(.x, 0)))

        method_output_dir <- file.path(OUTPUT_DIR, method)
        dir.create(method_output_dir, recursive = TRUE, showWarnings = FALSE)
        readr::write_csv(method_matrix, file.path(method_output_dir, "merged_matrix.csv"))
        cat(method, "->", nrow(method_matrix), "features x",
            ncol(method_matrix) - 1, "samples\n")
        rm(method_matrix)
    } else {
        rm(all_df)
        cat(method, ": no input files found\n")
    }
    gc()
}
