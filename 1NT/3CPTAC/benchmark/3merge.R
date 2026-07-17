rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)
library(tidyr)

# 1. Basic settings
BASE_DIR <- paste0(PROJ_ROOT, "/1NT/3CPTAC/benchmark")
LINKS_FILE <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")


available_datasets <- c(
"LUAD", "eLUAD", "KIRC", "PRAD"
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

# SSN
method <- "SSN"
all_df <- list()

for (dataset in available_datasets) {
    file <- file.path(BASE_DIR, method, dataset, "delta.csv")
    if (!file.exists(file)) next

    one_dataset <- readr::read_csv(file, show_col_types = FALSE)
    names(one_dataset)[1:2] <- c("protein1", "protein2")

    one_dataset <- one_dataset %>%
        tidyr::pivot_longer(
            cols = -c(protein1, protein2),
            names_to = "Sample",
            values_to = "Weight"
        )

    all_df[[dataset]] <- one_dataset
}

ssn_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(ssn_long) > 0) {
    ssn_matrix <- make_matrix_from_long(ssn_long)
    rm(ssn_long)
    readr::write_csv(ssn_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(ssn_matrix), "interactions x", ncol(ssn_matrix) - 1, "samples\n")
    rm(ssn_matrix)
} else {
    rm(ssn_long)
    cat(method, ": no input files found\n")
}
gc()

# ===================================================

# ===================================================


# SWEET
method <- "SWEET"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, method, dataset)
    files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    files <- files[!grepl("zscore|mean_std|weight|runtime", basename(files))]
    if (length(files) == 0) next

    one_dataset <- lapply(files, function(f) {
        sample_id <- sub("\\.txt$", "", basename(f))
        read.table(f, header = TRUE, fill = TRUE) %>%
            dplyr::transmute(
                protein1 = gene1,
                protein2 = gene2,
                Sample = sample_id,
                Weight = raw_edge_score
            )
    }) %>%
        dplyr::bind_rows()

    all_df[[dataset]] <- one_dataset
}

sweet_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(sweet_long) > 0) {
    sweet_matrix <- make_matrix_from_long(sweet_long)
    rm(sweet_long)
    readr::write_csv(sweet_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(sweet_matrix), "interactions x", ncol(sweet_matrix) - 1, "samples\n")
    rm(sweet_matrix)
} else {
    rm(sweet_long)
    cat(method, ": no input files found\n")
}
gc()


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

mossn_uniform_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(mossn_uniform_long) > 0) {
    mossn_uniform_matrix <- make_matrix_from_long(mossn_uniform_long)
    rm(mossn_uniform_long)
    readr::write_csv(mossn_uniform_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(mossn_uniform_matrix), "interactions x", ncol(mossn_uniform_matrix) - 1, "samples\n")
    rm(mossn_uniform_matrix)
} else {
    rm(mossn_uniform_long)
    cat(method, ": no input files found\n")
}
gc()

# ===================================================

# ===================================================


# Patkar
method <- "Patkar"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, method, dataset)
    files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    if (length(files) == 0) next

    one_dataset <- lapply(files, function(f) {
        sample_id <- sub("\\.txt$", "", basename(f))
        readr::read_tsv(f, show_col_types = FALSE) %>%
            dplyr::transmute(
                protein1 = gene1,
                protein2 = gene2,
                Sample = sample_id,
                Weight = score
            )
    }) %>%
        dplyr::bind_rows()

    all_df[[dataset]] <- one_dataset
}

patkar_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(patkar_long) > 0) {
    patkar_matrix <- make_matrix_from_long(patkar_long)
    rm(patkar_long)
    readr::write_csv(patkar_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(patkar_matrix), "interactions x", ncol(patkar_matrix) - 1, "samples\n")
    rm(patkar_matrix)
} else {
    rm(patkar_long)
    cat(method, ": no input files found\n")
}
gc()

# ===================================================

# ===================================================


# PPIXpress
method <- "PPIXpress"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, method, dataset)
    files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    if (length(files) == 0) next

    one_dataset <- lapply(files, function(f) {
        sample_id <- sub("\\.txt$", "", basename(f))
        readr::read_tsv(f, show_col_types = FALSE) %>%
            dplyr::transmute(
                protein1 = gene1,
                protein2 = gene2,
                Sample = sample_id,
                Weight = score
            )
    }) %>%
        dplyr::bind_rows()

    all_df[[dataset]] <- one_dataset
}

ppixpress_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(ppixpress_long) > 0) {
    ppixpress_matrix <- make_matrix_from_long(ppixpress_long)
    rm(ppixpress_long)
    readr::write_csv(ppixpress_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(ppixpress_matrix), "interactions x", ncol(ppixpress_matrix) - 1, "samples\n")
    rm(ppixpress_matrix)
} else {
    rm(ppixpress_long)
    cat(method, ": no input files found\n")
}
gc()

# ===================================================

# ===================================================


# Proteinarium
method <- "Proteinarium"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, method, dataset)
    files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
    if (length(files) == 0) next

    one_dataset <- lapply(files, function(f) {
        sample_id <- sub("\\.txt$", "", basename(f))
        readr::read_tsv(f, show_col_types = FALSE) %>%
            dplyr::transmute(
                protein1 = gene1,
                protein2 = gene2,
                Sample = sample_id,
                Weight = score
            )
    }) %>%
        dplyr::bind_rows()

    all_df[[dataset]] <- one_dataset
}

proteinarium_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(proteinarium_long) > 0) {
    proteinarium_matrix <- make_matrix_from_long(proteinarium_long)
    rm(proteinarium_long)
    readr::write_csv(proteinarium_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(proteinarium_matrix), "interactions x", ncol(proteinarium_matrix) - 1, "samples\n")
    rm(proteinarium_matrix)
} else {
    rm(proteinarium_long)
    cat(method, ": no input files found\n")
}
gc()

# ===================================================

# ===================================================


# LIONESS
method <- "LIONESS"
all_df <- list()

for (dataset in available_datasets) {
    file <- file.path(BASE_DIR, method, dataset, "result.csv")
    if (!file.exists(file)) next

    one_dataset <- readr::read_csv(file, show_col_types = FALSE)
    names(one_dataset)[1:2] <- c("protein1", "protein2")

    one_dataset <- one_dataset %>%
        tidyr::pivot_longer(
            cols = -c(protein1, protein2),
            names_to = "Sample",
            values_to = "Weight"
        )

    all_df[[dataset]] <- one_dataset
}

lioness_long <- dplyr::bind_rows(all_df)
rm(all_df)
if (nrow(lioness_long) > 0) {
    lioness_matrix <- make_matrix_from_long(lioness_long)
    rm(lioness_long)
    readr::write_csv(lioness_matrix, file.path(BASE_DIR, method, "merged_matrix.csv"))
    cat(method, "->", nrow(lioness_matrix), "interactions x", ncol(lioness_matrix) - 1, "samples\n")
    rm(lioness_matrix)
} else {
    rm(lioness_long)
    cat(method, ": no input files found\n")
}
gc()
