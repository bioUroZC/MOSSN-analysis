rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)

RESULTS_DIR <- paste0(PROJ_ROOT, "/4coupled/results")
CANCERS     <- c("ACC", "BLCA", "BRCA", "CESC", 
                 "CRC", "ESCA", "GBM", "HNSC", "KIRC",
                 "LGG", "LIHC", "LUAD", "LUSC", 
                 "PAAD", "PRAD", "SARC", "STAD")

METHOD_DIRS <- list(
    "MOSSN_EXP" = c("MOSSN_EXP", "EXP_single"),
    "MOSSN_MET" = c("MOSSN_MET", "MET_single"),
    "MOSSN_CNV" = c("MOSSN_CNV", "CNV_single"),
    "MOSSN_NoCross" = c("MOSSN_NoCross", "MUL_noCross"),
    "MOSSN_Restart" = c("MOSSN_Restart", "MUL_full"),
    "MOSSN_Direct" = c("MOSSN_Direct", "MUL_direct"),
    "MOSSN_DirectNoDyn" = c("MOSSN_DirectNoDyn", "MUL_direct_fixed"),
    "MOSSN_MultiLayer" = c("MOSSN_MultiLayer", "MUL_multilayer")
)

for (method in names(METHOD_DIRS)) {
    cat("\n[", method, "]\n", sep = "")
    cancer_matrices <- list()

    for (cancer in CANCERS) {
        folders <- file.path(RESULTS_DIR, METHOD_DIRS[[method]], cancer)
        folder <- folders[file.exists(folders)][1]
        if (is.na(folder)) {
            cat("  ", cancer, ": no folder\n")
            next
        }
        files  <- list.files(folder, pattern = "_edges\\.csv$", full.names = TRUE)
        if (length(files) == 0) { cat("  ", cancer, ": no files\n"); next }

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
        rm(ref, mat); gc()
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
        rm(cmat); gc()
    }

    dir.create(file.path(RESULTS_DIR, method), showWarnings = FALSE, recursive = TRUE)
    out_df   <- tibble::rownames_to_column(as.data.frame(full_mat, check.names = FALSE),
                                           var = "Interaction")
    out_file <- file.path(RESULTS_DIR, method, "merged_matrix.csv")
    write_csv(out_df, out_file)
    cat("  Saved ->", out_file, "\n")
    rm(full_mat, out_df); gc()
}
