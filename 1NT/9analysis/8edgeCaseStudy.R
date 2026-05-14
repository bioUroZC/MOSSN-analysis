rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
module_dir <- file.path(base_dir, "module_results")
out_root <- file.path(base_dir, "case_survival")

args <- commandArgs(trailingOnly = TRUE)
edge_name <- "MMP9_SPP1"

edge_arg <- args[grepl("^--edge=", args)]
if (length(edge_arg) > 0) {
  edge_name <- sub("^--edge=", "", edge_arg[1])
}

edge_dir_name <- gsub("[^A-Za-z0-9_]+", "_", edge_name)
out_dir <- file.path(out_root, edge_dir_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

resolve_edge_name <- function(edge_label, edge_pool) {
  if (edge_label %in% edge_pool) {
    return(edge_label)
  }
  parts <- strsplit(edge_label, "_", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    rev_label <- paste(rev(parts), collapse = "_")
    if (rev_label %in% edge_pool) {
      return(rev_label)
    }
  }
  NA_character_
}

recurrent_tbl <- read.csv(
  file.path(base_dir, "universal_recurrent_links.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

module_edges <- read.csv(
  file.path(module_dir, "module_edges.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
annotation_tbl <- read.csv(
  file.path(module_dir, "module_annotation_table.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

matched_edge <- resolve_edge_name(edge_name, unique(module_edges$link))
if (is.na(matched_edge)) {
  matched_edge <- resolve_edge_name(edge_name, unique(recurrent_tbl$link))
}
if (is.na(matched_edge)) {
  stop("Selected edge was not found in module or recurrence outputs: ", edge_name)
}

recurrent_row <- recurrent_tbl %>%
  dplyr::filter(link == matched_edge)

module_row <- module_edges %>%
  dplyr::filter(link == matched_edge) %>%
  dplyr::slice(1)

if (nrow(module_row) == 0) {
  stop("Selected edge was not assigned to a module: ", matched_edge)
}

annot_row <- annotation_tbl %>%
  dplyr::filter(module_id == module_row$module_id) %>%
  dplyr::slice(1)

module_context <- module_edges %>%
  dplyr::filter(module_id == module_row$module_id) %>%
  dplyr::arrange(dplyr::desc(weight), dplyr::desc(recurrent_count), link) %>%
  dplyr::mutate(is_selected_edge = link == matched_edge)

summary_tbl <- dplyr::bind_cols(
  data.frame(
    edge_query = edge_name,
    edge_matched = matched_edge,
    stringsAsFactors = FALSE
  ),
  recurrent_row %>%
    dplyr::select(
      link,
      n_gain,
      n_loss,
      n_cancers,
      recurrent_count,
      dominant_direction,
      consistency,
      recurrent_class,
      cancers_gain,
      cancers_loss,
      median_delta_across_cancers
    ),
  module_row %>%
    dplyr::select(module_id, direction, module_rank, weight),
  annot_row %>%
    dplyr::select(
      best_source_hk,
      best_term_hk,
      concise_label,
      hallmark_top,
      kegg_top,
      top_hubs
    )
)

utils::write.csv(
  summary_tbl,
  file.path(out_dir, "edge_case_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  module_context,
  file.path(out_dir, "edge_case_module_edges.csv"),
  row.names = FALSE
)

plot_df <- module_context %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::mutate(
    edge_label = factor(link, levels = rev(link)),
    fill_group = dplyr::if_else(is_selected_edge, "selected", "context")
  )

subtitle_txt <- paste0(
  module_row$module_id[1], " | ",
  annot_row$best_source_hk[1], " | ",
  annot_row$concise_label[1]
)

p <- ggplot(plot_df, aes(x = weight, y = edge_label, fill = fill_group)) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c(selected = "#B2182B", context = "#bdbdbd"),
    breaks = c("selected", "context"),
    labels = c("Selected edge", "Module context")
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Module edge weight",
    y = "Edges in selected module",
    fill = NULL,
    title = paste0("Case-study edge context: ", matched_edge),
    subtitle = subtitle_txt
  )

pdf(file.path(out_dir, "edge_case_module_context.pdf"), width = 10, height = 6)
print(p)
dev.off()

message("Saved edge case-study outputs to: ", out_dir)
