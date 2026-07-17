rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(data.table)
library(dplyr)

base_dir <- paste0(PROJ_ROOT, "/1NT/7casestudy")
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- paste0(PROJ_ROOT, "/1NT/2matrix/MOSS_full/merged_matrix.csv")
expr_file <- paste0(PROJ_ROOT, "/1NT/5parameter/data/LUAD_paired_expr.csv")
metadata_file <- paste0(PROJ_ROOT, "/1NT/5parameter/data/LUAD_paired_metadata.csv")
link_file <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")

top_edge_n <- 30
top_node_label_n <- 10

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
links <- read.csv(link_file, row.names = 1, stringsAsFactors = FALSE)
pair_matrix <- fread(matrix_file) |> as.data.frame()

if ("V1" %in% colnames(pair_matrix)) pair_matrix$V1 <- NULL
if (!"Interaction" %in% colnames(pair_matrix)) {
  names(pair_matrix)[1] <- "Interaction"
}

sample_ids <- metadata$sample
pair_matrix <- pair_matrix[, c("Interaction", intersect(sample_ids, colnames(pair_matrix))), drop = FALSE]

patient_ids <- unique(metadata$patient)

patient_distance <- lapply(patient_ids, function(pid) {
  tumor_sample <- metadata$sample[metadata$patient == pid & metadata$group == "01A"]
  normal_sample <- metadata$sample[metadata$patient == pid & metadata$group == "11A"]

  if (length(tumor_sample) == 0 || length(normal_sample) == 0) {
    return(NULL)
  }

  tumor_values <- pair_matrix[[tumor_sample[1]]]
  normal_values <- pair_matrix[[normal_sample[1]]]
  diff_values <- tumor_values - normal_values

  data.frame(
    patient = pid,
    tumor_sample = tumor_sample[1],
    normal_sample = normal_sample[1],
    edge_distance = sqrt(sum(diff_values^2, na.rm = TRUE)),
    mean_abs_edge_change = mean(abs(diff_values), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

patient_distance_df <- bind_rows(patient_distance)
write.csv(patient_distance_df, file.path(out_dir, "case_patient_distance.csv"), row.names = FALSE)

patient_distance_df <- patient_distance_df[order(patient_distance_df$edge_distance), , drop = FALSE]
mid_idx <- ceiling(nrow(patient_distance_df) / 2)
selected_case <- patient_distance_df[mid_idx, , drop = FALSE]

tumor_sample <- selected_case$tumor_sample[1]
normal_sample <- selected_case$normal_sample[1]
selected_patient <- selected_case$patient[1]

edge_change_df <- data.frame(
  Interaction = pair_matrix$Interaction,
  tumor_weight = pair_matrix[[tumor_sample]],
  normal_weight = pair_matrix[[normal_sample]],
  stringsAsFactors = FALSE
)
edge_change_df$delta <- edge_change_df$tumor_weight - edge_change_df$normal_weight
edge_change_df$abs_delta <- abs(edge_change_df$delta)
edge_change_df$direction <- ifelse(edge_change_df$delta > 0, "gain", "loss")

link_map <- links[, c("protein1", "protein2", "score")]
link_map$Interaction <- paste0(link_map$protein1, "_", link_map$protein2)
edge_change_df <- merge(edge_change_df, link_map, by = "Interaction", all.x = TRUE)
edge_change_df <- edge_change_df[order(edge_change_df$abs_delta, decreasing = TRUE), , drop = FALSE]

write.csv(edge_change_df, file.path(out_dir, "case_all_edge_changes.csv"), row.names = FALSE)

top_edge_df <- edge_change_df[seq_len(min(top_edge_n, nrow(edge_change_df))), , drop = FALSE]
write.csv(top_edge_df, file.path(out_dir, "case_top_edge_changes.csv"), row.names = FALSE)

display_nodes <- sort(unique(c(top_edge_df$protein1, top_edge_df$protein2)))
display_edge_df <- edge_change_df[edge_change_df$protein1 %in% display_nodes & edge_change_df$protein2 %in% display_nodes, , drop = FALSE]
display_edge_df$edge_width <- rank(display_edge_df$abs_delta, ties.method = "average") / max(1, nrow(display_edge_df))

node_stats <- bind_rows(
  data.frame(gene = display_edge_df$protein1, abs_delta = display_edge_df$abs_delta, direction = display_edge_df$direction, stringsAsFactors = FALSE),
  data.frame(gene = display_edge_df$protein2, abs_delta = display_edge_df$abs_delta, direction = display_edge_df$direction, stringsAsFactors = FALSE)
) |>
  group_by(gene) |>
  summarise(
    rewiring_degree = n(),
    rewiring_burden = sum(abs_delta, na.rm = TRUE),
    gain_edges = sum(direction == "gain", na.rm = TRUE),
    loss_edges = sum(direction == "loss", na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(rewiring_burden))

node_stats$label <- ifelse(seq_len(nrow(node_stats)) <= top_node_label_n, node_stats$gene, "")

write.csv(display_edge_df, file.path(out_dir, "case_display_edges.csv"), row.names = FALSE)
write.csv(node_stats, file.path(out_dir, "case_display_nodes.csv"), row.names = FALSE)

case_summary <- data.frame(
  patient = selected_patient,
  tumor_sample = tumor_sample,
  normal_sample = normal_sample,
  edge_distance_rank = mid_idx,
  total_patients = nrow(patient_distance_df),
  displayed_edges = nrow(display_edge_df),
  displayed_nodes = nrow(node_stats),
  top_changed_edges = nrow(top_edge_df),
  gain_edges_in_display = sum(display_edge_df$direction == "gain", na.rm = TRUE),
  loss_edges_in_display = sum(display_edge_df$direction == "loss", na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(case_summary, file.path(out_dir, "case_selection_summary.csv"), row.names = FALSE)

cat("Selected patient:", selected_patient, "\n")
cat("Tumor sample:", tumor_sample, "\n")
cat("Normal sample:", normal_sample, "\n")
cat("Saved outputs to:", out_dir, "\n")
