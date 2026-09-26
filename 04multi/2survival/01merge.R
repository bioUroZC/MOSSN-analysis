rm(list = ls())
library(dplyr)
library(readr)

RESULTS_DIR <- "/proj/c.zihao/work1/04multi/2survival/merge"
CANCERS     <- c("BLCA", "LIHC", "LUAD", "SARC", "STAD")

METHODS <- c(
    "EXP", "MET", "CNV",
    "NoCross", "Restart",
    "Direct", "NoDyn", "MultiLayer"
)

for (method in METHODS) {
    cat("\n[", method, "]\n", sep = "")
    cancer_matrices <- list()

    for (cancer in CANCERS) {
        folder <- file.path(RESULTS_DIR, method, cancer)
        if (!dir.exists(folder)) {
            cat("  ", cancer, ": no folder\n")
            next
        }
        files  <- list.files(folder, pattern = "_edges\\.csv$", full.names = TRUE)
        if (length(files) == 0) {
            cat("  ", cancer, ": no files\n")
            next
        }

        ref             <- read_csv(files[1], show_col_types = FALSE)
        interaction_ids <- paste(ref$Node1, ref$Node2, sep = "--")
        sample_ids      <- sub("_edges\\.csv$", "", basename(files))

        mat <- matrix(0, nrow = length(interaction_ids),
                         ncol = length(sample_ids),
                         dimnames = list(interaction_ids, sample_ids))

        for (i in seq_along(files)) {
            mat[, i] <- read_csv(files[i], show_col_types = FALSE,
                                 col_select = "FinalWeight")[[1]]
        }

        cancer_matrices[[cancer]] <- mat
        cat("  ", cancer, ":", length(interaction_ids), "interactions x",
            length(sample_ids), "samples\n")
        rm(ref, mat)
        gc()
    }

    if (length(cancer_matrices) == 0) next

    all_interactions <- sort(unique(unlist(lapply(cancer_matrices, rownames))))
    all_samples      <- unlist(lapply(cancer_matrices, colnames))
    cat("  Combining:", length(all_interactions), "interactions x",
        length(all_samples), "samples\n")

    full_mat <- matrix(0, nrow = length(all_interactions), ncol = length(all_samples),
                       dimnames = list(all_interactions, all_samples))

    col_ptr <- 1L
    for (cancer in names(cancer_matrices)) {
        cmat    <- cancer_matrices[[cancer]]
        row_idx <- match(rownames(cmat), all_interactions)
        col_end <- col_ptr + ncol(cmat) - 1L
        full_mat[row_idx, col_ptr:col_end] <- cmat
        col_ptr <- col_end + 1L
        cancer_matrices[[cancer]] <- NULL
        rm(cmat)
        gc()
    }

    dir.create(file.path(RESULTS_DIR, method), showWarnings = FALSE, recursive = TRUE)
    out_df   <- tibble::rownames_to_column(as.data.frame(full_mat, check.names = FALSE),
                                           var = "Interaction")
    out_file <- file.path(RESULTS_DIR, method, "merged_matrix.csv")
    write_csv(out_df, out_file)
    cat("  Saved ->", out_file, "\n")
    rm(full_mat, out_df)
    gc()
}
