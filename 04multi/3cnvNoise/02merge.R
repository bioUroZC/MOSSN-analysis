rm(list = ls())
library(dplyr)
library(readr)

CANCER       <- "LUAD"
RESULTS_DIR  <- "/proj/c.zihao/work1/04multi/3cnvNoise/03merge"
NOISE_LEVELS <- c("k0.0", "k0.5", "k1.0", "k1.5", "k2.0", "k3.0", "k5.0")
METHODS      <- c(
    "CNV", "Restart",
    "Direct", "NoDyn", "MultiLayer"
)

for (method in METHODS) {
    cat("\n[", method, "]\n", sep = "")

    for (klabel in NOISE_LEVELS) {
        folder <- file.path(RESULTS_DIR, method, klabel, CANCER)
        if (!dir.exists(folder)) {
            cat("  ", klabel, ": no folder, skip\n")
            next
        }
        files  <- list.files(folder, pattern = "_edges\\.csv$", full.names = TRUE)
        if (length(files) == 0) {
            cat("  ", klabel, ": no files, skip\n")
            next
        }

        ref             <- read_csv(files[1], show_col_types = FALSE)
        interaction_ids <- paste(ref$Node1, ref$Node2, sep = "--")
        sample_ids      <- sub("_edges\\.csv$", "", basename(files))

        mat <- matrix(0, nrow = length(interaction_ids), ncol = length(sample_ids),
                      dimnames = list(interaction_ids, sample_ids))

        for (i in seq_along(files)) {
            mat[, i] <- read_csv(files[i], show_col_types = FALSE,
                                 col_select = "FinalWeight")[[1]]
        }

        dir.create(file.path(RESULTS_DIR, method, klabel), showWarnings = FALSE, recursive = TRUE)
        out_df   <- tibble::rownames_to_column(as.data.frame(mat, check.names = FALSE),
                                               var = "Interaction")
        out_file <- file.path(RESULTS_DIR, method, klabel, "merged_matrix.csv")
        write_csv(out_df, out_file)
        cat("  ", klabel, ":", length(interaction_ids), "x", length(sample_ids),
            "-> merged_matrix.csv\n")
        rm(mat, out_df)
        gc()
    }
}
