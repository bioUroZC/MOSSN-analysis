rm(list = ls())

NOISE_DIR <- "/proj/c.zihao/work1/02analysis/02noise"
FULL_MATRIX_DIR <- "/proj/c.zihao/work1/01NT/01string/02benchmark/02merge/merged_matrices"
LEVELS <- c("05", "10", "15", "20")
METHODS <- c("MOSSN", "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress", "Proteinarium")

load_matrix <- function(path) {
  as.matrix(read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
}

for (level in LEVELS) {
  level_dir <- file.path(NOISE_DIR, level)
  selected_file <- file.path(level_dir, "data", "LUAD_noise_samples.csv")
  if (!file.exists(selected_file)) stop(sprintf("Missing selected-sample file: %s", selected_file))
  selected_samples <- read.csv(selected_file, stringsAsFactors = FALSE)$sample

  summary_rows <- list()
  sample_rows <- list()

  for (method in METHODS) {
    noisy_file <- file.path(level_dir, "merged_matrices", method, "merged_matrix.csv")
    full_file <- file.path(FULL_MATRIX_DIR, method, "merged_matrix.csv")
    if (!file.exists(noisy_file) || !file.exists(full_file)) {
      warning(sprintf("%s (%s%% noise): missing noisy or full matrix", method, level))
      next
    }

    noisy_matrix <- load_matrix(noisy_file)
    full_matrix <- load_matrix(full_file)
    common_edges <- intersect(rownames(noisy_matrix), rownames(full_matrix))
    common_samples <- selected_samples[
      selected_samples %in% colnames(noisy_matrix) & selected_samples %in% colnames(full_matrix)
    ]
    if (length(common_edges) == 0 || length(common_samples) == 0) {
      warning(sprintf("%s (%s%% noise): no common edges or samples", method, level))
      next
    }

    noisy_matrix <- noisy_matrix[common_edges, common_samples, drop = FALSE]
    full_matrix <- full_matrix[common_edges, common_samples, drop = FALSE]
    sample_spearman <- vapply(common_samples, function(sample_id) {
      suppressWarnings(cor(
        as.numeric(full_matrix[, sample_id]),
        as.numeric(noisy_matrix[, sample_id]),
        method = "spearman", use = "pairwise.complete.obs"
      ))
    }, numeric(1))

    sample_rows[[length(sample_rows) + 1]] <- data.frame(
      noise_pct = as.integer(level), method = method, sample = common_samples,
      spearman = sample_spearman, stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      noise_pct = as.integer(level), method = method,
      consistency = mean(sample_spearman, na.rm = TRUE),
      sd_spearman = stats::sd(sample_spearman, na.rm = TRUE),
      n_samples = length(common_samples), n_edges = length(common_edges),
      stringsAsFactors = FALSE
    )
  }

  output_dir <- file.path(level_dir, "consistency")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
  sample_df <- if (length(sample_rows) > 0) do.call(rbind, sample_rows) else data.frame()
  if (nrow(summary_df) > 0) {
    summary_df <- summary_df[order(summary_df$consistency, decreasing = TRUE), , drop = FALSE]
  }
  write.csv(summary_df, file.path(output_dir, "consistency_df.csv"), row.names = FALSE)
  write.csv(sample_df, file.path(output_dir, "sample_level_consistency.csv"), row.names = FALSE)
  cat(sprintf("%s%% noise: wrote consistency for %d methods to %s\n", level, nrow(summary_df), output_dir))
}
