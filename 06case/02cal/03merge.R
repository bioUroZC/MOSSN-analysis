rm(list = ls())

library(data.table)
library(dplyr)

base_dir <- "/proj/c.zihao/work1/06case"
edge_dir <- file.path(base_dir, "04results", "MOSSN")
out_file <- file.path(edge_dir, "merged_matrix.csv")
links_file <- "/proj/c.zihao/work1/01data/06string/links.csv"

edge_files <- list.files(edge_dir, pattern = "_edges\\.csv$", full.names = TRUE)
if (length(edge_files) == 0) {
  stop("No *_edges.csv files found in ", edge_dir)
}

links_ref <- fread(links_file) %>%
  dplyr::select(protein1, protein2) %>%
  dplyr::transmute(
    protein1 = pmin(protein1, protein2),
    protein2 = pmax(protein1, protein2)
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Interaction = paste(protein1, protein2, sep = "_"))

interaction_ids <- links_ref$Interaction
sample_ids <- sort(sub("_edges\\.csv$", "", basename(edge_files)))

out_mat <- matrix(
  0,
  nrow = length(interaction_ids),
  ncol = length(sample_ids),
  dimnames = list(interaction_ids, sample_ids)
)

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  one <- data.table::fread(
    file.path(edge_dir, paste0(sample_id, "_edges.csv")),
    select = c("Node1", "Node2", "FinalWeight")
  )
  one[, Interaction := paste(pmin(Node1, Node2), pmax(Node1, Node2), sep = "_")]
  row_id <- match(one$Interaction, interaction_ids)
  keep <- !is.na(row_id)
  out_mat[cbind(row_id[keep], i)] <- one$FinalWeight[keep]
  message(sample_id, ": ", sum(keep), " edges")
}

out_df <- data.frame(
  Interaction = interaction_ids,
  out_mat,
  check.names = FALSE,
  row.names = NULL
)
data.table::fwrite(out_df, out_file)

message(
  "Wrote ", out_file, ": ",
  nrow(out_df), " interactions x ",
  ncol(out_df) - 1, " samples"
)
