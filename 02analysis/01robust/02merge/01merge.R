rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)

ROBUST_DIR <- "/proj/c.zihao/work1/02analysis/01robust"
LEVELS <- c("10", "30", "50", "70")
DATASET <- "LUAD"
METHODS <- c("MOSSN", "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress", "Proteinarium")
LINKS_FILE <- "/proj/c.zihao/work1/01data/06string/links.csv"

links_ref <- readr::read_csv(LINKS_FILE, show_col_types = FALSE) %>%
  transmute(
    protein1_raw = protein1,
    protein2_raw = protein2,
    protein1 = pmin(protein1_raw, protein2_raw),
    protein2 = pmax(protein1_raw, protein2_raw)
  ) %>%
  select(protein1, protein2) %>%
  distinct()
links_ref$Interaction <- paste(links_ref$protein1, links_ref$protein2, sep = "_")

make_matrix_from_long <- function(df_long) {
  df_long <- df_long %>%
    transmute(
      protein1_raw = protein1,
      protein2_raw = protein2,
      protein1 = pmin(protein1_raw, protein2_raw),
      protein2 = pmax(protein1_raw, protein2_raw),
      Sample = Sample,
      Weight = as.numeric(Weight)
    )

  sample_ids <- sort(unique(df_long$Sample))
  interaction_ids <- links_ref$Interaction
  out_mat <- matrix(0, nrow = length(interaction_ids), ncol = length(sample_ids),
                    dimnames = list(interaction_ids, sample_ids))

  row_id <- match(paste(df_long$protein1, df_long$protein2, sep = "_"), interaction_ids)
  col_id <- match(df_long$Sample, sample_ids)
  keep <- !is.na(row_id) & !is.na(col_id) & !is.na(df_long$Weight)
  out_mat[cbind(row_id[keep], col_id[keep])] <- df_long$Weight[keep]

  out_df <- as.data.frame(out_mat, check.names = FALSE)
  out_df$Interaction <- rownames(out_mat)
  out_df[, c("Interaction", sample_ids), drop = FALSE]
}

read_method_results <- function(level_dir, method) {
  method_dir <- file.path(level_dir, method, DATASET)

  if (method == "MOSSN") {
    files <- list.files(method_dir, pattern = "_edges\\.csv$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    return(lapply(files, readr::read_csv, show_col_types = FALSE) %>%
      bind_rows() %>%
      transmute(protein1 = Node1, protein2 = Node2, Sample = Sample, Weight = FinalWeight))
  }

  if (method == "SSN") {
    file <- file.path(method_dir, "delta.csv")
    if (!file.exists(file)) return(NULL)
    wide <- readr::read_csv(file, show_col_types = FALSE)
    names(wide)[1:2] <- c("protein1", "protein2")
    return(pivot_longer(wide, cols = -c(protein1, protein2), names_to = "Sample", values_to = "Weight"))
  }

  if (method == "LIONESS") {
    file <- file.path(method_dir, "result.csv")
    if (!file.exists(file)) return(NULL)
    wide <- readr::read_csv(file, show_col_types = FALSE)
    names(wide)[1:2] <- c("protein1", "protein2")
    return(pivot_longer(wide, cols = -c(protein1, protein2), names_to = "Sample", values_to = "Weight"))
  }

  files <- list.files(method_dir, pattern = "\\.txt$", full.names = TRUE)
  if (method == "SWEET") {
    files <- files[!grepl("zscore|mean_std|weight|runtime", basename(files))]
  }
  if (length(files) == 0) return(NULL)

  if (method == "SWEET") {
    return(lapply(files, function(file) {
      read.table(file, header = TRUE, fill = TRUE) %>%
        transmute(
          protein1 = gene1,
          protein2 = gene2,
          Sample = sub("\\.txt$", "", basename(file)),
          Weight = raw_edge_score
        )
    }) %>% bind_rows())
  }

  lapply(files, function(file) {
    readr::read_tsv(file, show_col_types = FALSE) %>%
      transmute(
        protein1 = gene1,
        protein2 = gene2,
        Sample = sub("\\.txt$", "", basename(file)),
        Weight = score
      )
  }) %>% bind_rows()
}

for (level in LEVELS) {
  level_dir <- file.path(ROBUST_DIR, level)
  output_dir <- file.path(level_dir, "merged_matrices")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (method in METHODS) {
    long_df <- read_method_results(level_dir, method)
    if (is.null(long_df) || nrow(long_df) == 0) {
      warning(sprintf("%s (%s%%): no raw results found", method, level))
      next
    }

    method_dir <- file.path(output_dir, method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    matrix_df <- make_matrix_from_long(long_df)
    output_file <- file.path(method_dir, "merged_matrix.csv")
    readr::write_csv(matrix_df, output_file)
    cat(sprintf("%s (%s%%): %d interactions x %d samples -> %s\n",
                method, level, nrow(matrix_df), ncol(matrix_df) - 1, output_file))
  }
}
