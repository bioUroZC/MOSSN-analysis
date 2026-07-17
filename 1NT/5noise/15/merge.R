rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(readr)
library(tidyr)

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  match_idx <- grep("--file=", cmd_args)
  if (length(match_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[match_idx[1]]), mustWork = TRUE))
  }
  frame_files <- vapply(sys.frames(), function(env) {
    if (exists("ofile", envir = env, inherits = FALSE))
      get("ofile", envir = env, inherits = FALSE)
    else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0)
    return(normalizePath(frame_files[length(frame_files)], mustWork = TRUE))
  stop("Unable to determine script path.")
}

BASE_DIR    <- dirname(get_script_path())
DATASET     <- "LUAD"
LINKS_FILE  <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")

links_ref <- readr::read_csv(LINKS_FILE, show_col_types = FALSE) %>%
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
  df_long <- df_long %>%
    dplyr::transmute(
      protein1_raw = protein1,
      protein2_raw = protein2,
      protein1 = pmin(protein1_raw, protein2_raw),
      protein2 = pmax(protein1_raw, protein2_raw),
      Sample   = Sample,
      Weight   = Weight
    )
  sample_ids      <- sort(unique(df_long$Sample))
  interaction_ids <- links_ref$Interaction
  out_mat <- matrix(0, nrow = length(interaction_ids), ncol = length(sample_ids))
  rownames(out_mat) <- interaction_ids
  colnames(out_mat) <- sample_ids
  cur_int <- paste(df_long$protein1, df_long$protein2, sep = "_")
  row_id  <- match(cur_int, interaction_ids)
  col_id  <- match(df_long$Sample, sample_ids)
  keep    <- !is.na(row_id) & !is.na(col_id)
  out_mat[cbind(row_id[keep], col_id[keep])] <- df_long$Weight[keep]
  out_df <- as.data.frame(out_mat, check.names = FALSE)
  out_df$Interaction <- rownames(out_mat)
  out_df[, c("Interaction", sample_ids)]
}

save_merged <- function(method, df_long) {
  if (nrow(df_long) == 0) {
    cat(method, ": no input files found\n")
    return(invisible(NULL))
  }
  mat <- make_matrix_from_long(df_long)
  readr::write_csv(mat, file.path(BASE_DIR, method, "merged_matrix.csv"))
  cat(method, "->", nrow(mat), "interactions x", ncol(mat) - 1, "samples\n")
}

# ------------------------------------------------------------------
# MOSSN_uniform
method <- "MOSSN_uniform"
folder <- file.path(BASE_DIR, method, DATASET)
files <- list.files(folder, pattern = "_edges\\.csv$", full.names = TRUE)
long_df <- lapply(files, readr::read_csv, show_col_types = FALSE) %>%
  dplyr::bind_rows() %>%
  dplyr::transmute(protein1 = Node1, protein2 = Node2, Sample = Sample, Weight = FinalWeight)
save_merged(method, long_df)

# ------------------------------------------------------------------
# SSN
method  <- "SSN"
file    <- file.path(BASE_DIR, method, DATASET, "delta.csv")
if (file.exists(file)) {
  wide <- readr::read_csv(file, show_col_types = FALSE)
  names(wide)[1:2] <- c("protein1", "protein2")
  long_df <- tidyr::pivot_longer(wide, cols = -c(protein1, protein2),
                                 names_to = "Sample", values_to = "Weight")
} else {
  long_df <- dplyr::tibble()
}
save_merged(method, long_df)

# ------------------------------------------------------------------
# SWEET
method  <- "SWEET"
folder  <- file.path(BASE_DIR, method, DATASET)
files   <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
files   <- files[!grepl("zscore|mean_std|weight|runtime", basename(files))]
long_df <- lapply(files, function(f) {
  sample_id <- sub("\\.txt$", "", basename(f))
  read.table(f, header = TRUE, fill = TRUE) %>%
    dplyr::transmute(protein1 = gene1, protein2 = gene2,
                     Sample = sample_id, Weight = raw_edge_score)
}) %>% dplyr::bind_rows()
save_merged(method, long_df)

# ------------------------------------------------------------------
# LIONESS
method  <- "LIONESS"
file    <- file.path(BASE_DIR, method, DATASET, "result.csv")
if (file.exists(file)) {
  wide <- readr::read_csv(file, show_col_types = FALSE)
  names(wide)[1:2] <- c("protein1", "protein2")
  long_df <- tidyr::pivot_longer(wide, cols = -c(protein1, protein2),
                                 names_to = "Sample", values_to = "Weight")
} else {
  long_df <- dplyr::tibble()
}
save_merged(method, long_df)

# ------------------------------------------------------------------
# Patkar
method  <- "Patkar"
folder  <- file.path(BASE_DIR, method, DATASET)
files   <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
long_df <- lapply(files, function(f) {
  sample_id <- sub("\\.txt$", "", basename(f))
  readr::read_tsv(f, show_col_types = FALSE) %>%
    dplyr::transmute(protein1 = gene1, protein2 = gene2,
                     Sample = sample_id, Weight = score)
}) %>% dplyr::bind_rows()
save_merged(method, long_df)

# ------------------------------------------------------------------
# PPIXpress
method  <- "PPIXpress"
folder  <- file.path(BASE_DIR, method, DATASET)
files   <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
long_df <- lapply(files, function(f) {
  sample_id <- sub("\\.txt$", "", basename(f))
  readr::read_tsv(f, show_col_types = FALSE) %>%
    dplyr::transmute(protein1 = gene1, protein2 = gene2,
                     Sample = sample_id, Weight = score)
}) %>% dplyr::bind_rows()
save_merged(method, long_df)

# ------------------------------------------------------------------
# Proteinarium
method  <- "Proteinarium"
folder  <- file.path(BASE_DIR, method, DATASET)
files   <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
long_df <- lapply(files, function(f) {
  sample_id <- sub("\\.txt$", "", basename(f))
  readr::read_tsv(f, show_col_types = FALSE) %>%
    dplyr::transmute(protein1 = gene1, protein2 = gene2,
                     Sample = sample_id, Weight = score)
}) %>% dplyr::bind_rows()
save_merged(method, long_df)
