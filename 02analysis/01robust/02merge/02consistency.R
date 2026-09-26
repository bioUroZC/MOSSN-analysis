rm(list = ls())

ROBUST_DIR <- "/proj/c.zihao/work1/02analysis/01robust"
FULL_MATRIX_DIR <- "/proj/c.zihao/work1/01NT/01string/02benchmark/02merge/merged_matrices"
LEVELS <- c("10", "30", "50", "70")
METHODS <- c("MOSSN", "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress", "Proteinarium")

load_matrix <- function(path) {
  matrix_df <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  as.matrix(matrix_df)
}

for (level in LEVELS) {
  level_dir <- file.path(ROBUST_DIR, level)
  selected_file <- file.path(level_dir, "data", "LUAD_half_samples.csv")
  if (!file.exists(selected_file)) {
    stop(sprintf("Missing selected-sample file: %s", selected_file))
  }
  selected_samples <- read.csv(selected_file, stringsAsFactors = FALSE)$sample

  summary_rows <- list()
  sample_rows <- list()

  for (method in METHODS) {
    subset_file <- file.path(level_dir, "merged_matrices", method, "merged_matrix.csv")
    full_file <- file.path(FULL_MATRIX_DIR, method, "merged_matrix.csv")
    if (!file.exists(subset_file) || !file.exists(full_file)) {
      warning(sprintf("%s (%s%%): missing subset or full matrix", method, level))
      next
    }

    subset_matrix <- load_matrix(subset_file)
    full_matrix <- load_matrix(full_file)
    common_edges <- intersect(rownames(subset_matrix), rownames(full_matrix))
    common_samples <- selected_samples[
      selected_samples %in% colnames(subset_matrix) & selected_samples %in% colnames(full_matrix)
    ]
    if (length(common_edges) == 0 || length(common_samples) == 0) {
      warning(sprintf("%s (%s%%): no common edges or samples", method, level))
      next
    }

    subset_matrix <- subset_matrix[common_edges, common_samples, drop = FALSE]
    full_matrix <- full_matrix[common_edges, common_samples, drop = FALSE]
    sample_spearman <- vapply(common_samples, function(sample_id) {
      suppressWarnings(cor(
        as.numeric(full_matrix[, sample_id]),
        as.numeric(subset_matrix[, sample_id]),
        method = "spearman", use = "pairwise.complete.obs"
      ))
    }, numeric(1))

    sample_rows[[length(sample_rows) + 1]] <- data.frame(
      retention_pct = as.integer(level), method = method, sample = common_samples,
      spearman = sample_spearman, stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      retention_pct = as.integer(level), method = method,
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
  cat(sprintf("%s%%: wrote consistency for %d methods to %s\n", level, nrow(summary_df), output_dir))
}
