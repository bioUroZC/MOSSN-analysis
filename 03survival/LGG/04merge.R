rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)

BASE_DIR <- "/proj/c.zihao/work1/03survival/LGG/"
OUTPUT_DIR <- "/proj/c.zihao/work1/03survival/LGG/02matrix"
LINKS_FILE <- "/proj/c.zihao/work1/01data/06string/links.csv"

available_datasets <- c("CGGA301", "CGGA325", "GSE16011", "TCGALGG")

EDGE_METHODS <- c("MOSSN", "noCorr", "noRWR", "RandomBackbone")
FEATURE_METHODS <- c("RawExpr", "NodeRWR")

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

reference_genes <- sort(unique(c(links_ref$protein1, links_ref$protein2)))

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

write_interaction_matrix <- function(method, df_long) {
    if (nrow(df_long) == 0) {
        cat(method, ": no input files found\n")
        return(invisible(NULL))
    }
    mat <- make_matrix_from_long(df_long)
    readr::write_csv(mat, file.path(OUTPUT_DIR, paste0(method, ".csv")))
    cat(method, "->", nrow(mat), "interactions x", ncol(mat) - 1, "samples\n")
    invisible(NULL)
}

method <- "PPIXpress"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, dataset, method)
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
write_interaction_matrix(method, ppixpress_long)
rm(ppixpress_long)
gc()

method <- "Patkar"
all_df <- list()

for (dataset in available_datasets) {
    folder <- file.path(BASE_DIR, dataset, method)
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
write_interaction_matrix(method, patkar_long)
rm(patkar_long)
gc()

for (method in EDGE_METHODS) {
    all_df <- list()

    for (dataset in available_datasets) {
        folder <- file.path(BASE_DIR, dataset, method)
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

    edge_long <- dplyr::bind_rows(all_df)
    rm(all_df)
    write_interaction_matrix(method, edge_long)
    rm(edge_long)
    gc()
}

for (method in FEATURE_METHODS) {
    all_df <- list()

    for (dataset in available_datasets) {
        file <- file.path(BASE_DIR, dataset, method, paste0(dataset, ".csv"))
        if (!file.exists(file)) next

        one_dataset <- readr::read_csv(
            file,
            show_col_types = FALSE,
            name_repair = "minimal"
        )
        names(one_dataset)[1] <- "Interaction"
        one_dataset$Interaction <- as.character(one_dataset$Interaction)
        all_df[[dataset]] <- one_dataset
    }

    if (length(all_df) > 0) {
        merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Interaction"), all_df)
        ref <- data.frame(Interaction = reference_genes, stringsAsFactors = FALSE)
        feature_matrix <- dplyr::left_join(ref, merged, by = "Interaction")
        feature_matrix[is.na(feature_matrix)] <- 0
        sample_ids <- sort(setdiff(colnames(feature_matrix), "Interaction"))
        feature_matrix <- feature_matrix[, c("Interaction", sample_ids)]
        readr::write_csv(feature_matrix, file.path(OUTPUT_DIR, paste0(method, ".csv")))
        cat(method, "->", nrow(feature_matrix), "features x", ncol(feature_matrix) - 1, "samples\n")
        rm(merged, ref, feature_matrix)
    } else {
        cat(method, ": no input files found\n")
    }
    rm(all_df)
    gc()
}
