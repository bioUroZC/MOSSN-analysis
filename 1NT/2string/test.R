rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)

# 1. Basic settings
BASE_DIR <- "/proj/c.zihao/work1/1NT/2string"
LINKS_FILE <- "/proj/c.zihao/work1/1NT/1data/string/links.csv"


available_datasets <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

# 2. Read reference interactions
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

# 3. A small helper: convert long table to matrix
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

# ===================================================

# ===================================================


# MOSSN_noPrior
method <- "MOSSN_noPrior"
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

moss_noPrior_long <- dplyr::bind_rows(all_df)
if (nrow(moss_noPrior_long) > 0) {
    moss_noPrior_matrix <- make_matrix_from_long(moss_noPrior_long)
    readr::write_csv(moss_noPrior_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(moss_noPrior_matrix), "interactions x", ncol(moss_noPrior_matrix) - 1, "samples\n")
} else {
    cat(method, ": no input files found\n")
}

# ===================================================

# ===================================================


# MOSSN_noSeed
method <- "MOSSN_noSeed"
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

moss_noSeed_long <- dplyr::bind_rows(all_df)
if (nrow(moss_noSeed_long) > 0) {
    moss_noSeed_matrix <- make_matrix_from_long(moss_noSeed_long)
    readr::write_csv(moss_noSeed_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(moss_noSeed_matrix), "interactions x", ncol(moss_noSeed_matrix) - 1, "samples\n")
} else {
    cat(method, ": no input files found\n")
}



# ===================================================

# ===================================================

# MOSSN_noCorr
method <- "MOSSN_noCorr"
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

moss_noCorr_long <- dplyr::bind_rows(all_df)

if (nrow(moss_noCorr_long) > 0) {
  moss_noCorr_matrix <- make_matrix_from_long(moss_noCorr_long)
  readr::write_csv(moss_noCorr_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
  cat(method, "->", nrow(moss_noCorr_matrix), "interactions x", ncol(moss_noCorr_matrix) - 1, "samples\n")
} else {
  cat(method, ": no input files found\n")
}

# ===================================================

# ===================================================


# MOSSN_noRWR
method <- "MOSSN_noRWR"
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

moss_noRWR_long <- dplyr::bind_rows(all_df)
if (nrow(moss_noRWR_long) > 0) {
    moss_noRWR_matrix <- make_matrix_from_long(moss_noRWR_long)
    readr::write_csv(moss_noRWR_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(moss_noRWR_matrix), "interactions x", ncol(moss_noRWR_matrix) - 1, "samples\n")
} else {
    cat(method, ": no input files found\n")
}

# ===================================================

# ===================================================


# MOSSN_uniform
method <- "MOSSN_uniform"
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

moss_uniform_long <- dplyr::bind_rows(all_df)
if (nrow(moss_uniform_long) > 0) {
    moss_uniform_matrix <- make_matrix_from_long(moss_uniform_long)
    readr::write_csv(moss_uniform_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(moss_uniform_matrix), "interactions x", ncol(moss_uniform_matrix) - 1, "samples\n")
} else {
    cat(method, ": no input files found\n")
}
