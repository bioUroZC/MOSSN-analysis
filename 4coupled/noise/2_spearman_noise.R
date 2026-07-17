rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)

CANCER       <- "LUAD"
RESULTS_DIR  <- paste0(PROJ_ROOT, "/4coupled/noise/results")
OUT_DIR      <- file.path(RESULTS_DIR, "spearman_noise")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NOISE_LEVELS <- c("k0.5", "k1.0", "k1.5", "k2.0", "k3.0", "k5.0")   # compare against k0.0
METHOD_DIRS  <- list(
    "MOSSN_EXP" = c("MOSSN_EXP", "EXP_single"),
    "MOSSN_NoCross" = c("MOSSN_NoCross", "MUL_noCross"),
    "MOSSN_Restart" = c("MOSSN_Restart", "MUL_full"),
    "MOSSN_Direct" = c("MOSSN_Direct", "MUL_direct"),
    "MOSSN_DirectNoDyn" = c("MOSSN_DirectNoDyn", "MOSSN_DirectFixed", "MUL_direct_fixed", "MUL_direct_nodyn"),
    "MOSSN_MultiLayer" = c("MOSSN_MultiLayer", "MUL_multilayer")
)
METHODS <- names(METHOD_DIRS)

results <- data.frame()

for (method in METHODS) {
    cat(sprintf("\n[%s]\n", method))

    ref_files <- file.path(RESULTS_DIR, METHOD_DIRS[[method]], "k0.0", "merged_matrix.csv")
    ref_file <- ref_files[file.exists(ref_files)][1]
    if (is.na(ref_file)) {
        cat("  k0.0 merged_matrix.csv not found, skip\n")
        next
    }
    ref <- read_csv(ref_file, show_col_types = FALSE) |> as.data.frame()
    rownames(ref) <- ref$Interaction
    ref$Interaction <- NULL

    for (klabel in NOISE_LEVELS) {
        noisy_files <- file.path(RESULTS_DIR, METHOD_DIRS[[method]], klabel, "merged_matrix.csv")
        noisy_file <- noisy_files[file.exists(noisy_files)][1]
        if (is.na(noisy_file)) {
            cat("  ", klabel, ": not found, skip\n")
            next
        }
        noisy <- read_csv(noisy_file, show_col_types = FALSE) |> as.data.frame()
        rownames(noisy) <- noisy$Interaction
        noisy$Interaction <- NULL

        common_samples      <- intersect(colnames(ref), colnames(noisy))
        common_interactions <- intersect(rownames(ref), rownames(noisy))

        if (length(common_samples) < 5 || length(common_interactions) < 10) {
            cat("  ", klabel, ": too few overlap, skip\n")
            next
        }

        ref_sub   <- ref[common_interactions,   common_samples]
        noisy_sub <- noisy[common_interactions, common_samples]

        # per-sample Spearman between original and noisy edge weight vectors
        rho_vec <- sapply(common_samples, function(s) {
            cor(ref_sub[[s]], noisy_sub[[s]], method = "spearman", use = "complete.obs")
        })

        results <- rbind(results, data.frame(
            Method        = method,
            NoiseLevel    = klabel,
            N_samples     = length(common_samples),
            N_interactions= length(common_interactions),
            Mean_Spearman = round(mean(rho_vec, na.rm = TRUE), 4),
            Median_Spearman = round(median(rho_vec, na.rm = TRUE), 4),
            SD_Spearman   = round(sd(rho_vec, na.rm = TRUE), 4),
            Q25_Spearman  = round(quantile(rho_vec, 0.25, na.rm = TRUE), 4),
            Q75_Spearman  = round(quantile(rho_vec, 0.75, na.rm = TRUE), 4)
        ))

        cat(sprintf("  %s | n=%d | median rho=%.4f (%.4f, %.4f)\n",
                    klabel, length(common_samples),
                    median(rho_vec, na.rm = TRUE),
                    quantile(rho_vec, 0.25, na.rm = TRUE),
                    quantile(rho_vec, 0.75, na.rm = TRUE)))
    }
}

write.csv(results, file.path(OUT_DIR, "spearman_results.csv"), row.names = FALSE)
cat("\nSaved ->", file.path(OUT_DIR, "spearman_results.csv"), "\n")
