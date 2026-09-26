rm(list = ls())
library(dplyr)
library(readr)

CANCER       <- "LUAD"
RESULTS_DIR  <- "/proj/c.zihao/work1/04multi/3metNoise/03merge"
OUT_DIR      <- "/proj/c.zihao/work1/04multi/3metNoise/04result"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NOISE_LEVELS <- c("k0.5", "k1.0", "k1.5", "k2.0", "k3.0", "k5.0")   
METHODS      <- c(
    "MET", "Restart",
    "Direct", "NoDyn", "MultiLayer"
)

results <- data.frame()

for (method in METHODS) {
    cat(sprintf("\n[%s]\n", method))

    ref_file <- file.path(RESULTS_DIR, method, "k0.0", "merged_matrix.csv")
    if (!file.exists(ref_file)) {
        cat("  k0.0 merged_matrix.csv not found, skip\n")
        next
    }
    ref <- read_csv(ref_file, show_col_types = FALSE) |> as.data.frame()
    rownames(ref) <- ref$Interaction
    ref$Interaction <- NULL

    for (klabel in NOISE_LEVELS) {
        noisy_file <- file.path(RESULTS_DIR, method, klabel, "merged_matrix.csv")
        if (!file.exists(noisy_file)) {
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

write.csv(results, file.path(OUT_DIR, "spearman.csv"), row.names = FALSE)
cat("\nSaved ->", file.path(OUT_DIR, "spearman.csv"), "\n")
