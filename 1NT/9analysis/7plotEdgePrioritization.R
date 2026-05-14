rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
module_dir <- file.path(base_dir, "module_results")

edge_file <- file.path(module_dir, "module_edges.csv")
annot_file <- file.path(module_dir, "module_annotation_table.csv")
out_csv <- file.path(base_dir, "edge_prioritization_candidates.csv")
out_pdf <- file.path(base_dir, "edge_prioritization_scatter.pdf")

target_edge <- "MMP9_SPP1"
label_n <- 6

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

rescale_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) < 1e-8) {
    return(rep(1, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

edges <- read.csv(
  edge_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ann <- read.csv(
  annot_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

matched_edge <- resolve_edge_name(target_edge, unique(edges$link))

plot_tbl <- edges %>%
  dplyr::left_join(
    ann %>%
      dplyr::select(module_id, best_source_hk, concise_label),
    by = "module_id"
  ) %>%
  dplyr::mutate(
    module_label = ifelse(
      is.na(concise_label) | concise_label == "",
      module_id,
      paste0(module_id, " | ", concise_label)
    ),
    priority_score = 0.5 * rescale_01(recurrent_count) +
      0.5 * rescale_01(weight),
    is_target = link == matched_edge,
    is_annotated = !is.na(concise_label) & concise_label != "" & concise_label != "Unannotated"
  ) %>%
  dplyr::arrange(dplyr::desc(priority_score), dplyr::desc(recurrent_count), dplyr::desc(weight), link)

utils::write.csv(plot_tbl, out_csv, row.names = FALSE)

label_tbl <- dplyr::bind_rows(
  plot_tbl %>% dplyr::filter(is_target),
  plot_tbl %>%
    dplyr::filter(is_annotated) %>%
    dplyr::slice_head(n = label_n)
) %>%
  dplyr::distinct(link, .keep_all = TRUE) %>%
  dplyr::mutate(label = link)

top_tbl <- plot_tbl %>%
  dplyr::filter(link %in% label_tbl$link)

plot_tbl <- plot_tbl %>%
  dplyr::mutate(
    highlight_group = dplyr::case_when(
      is_target ~ "Target edge",
      link %in% top_tbl$link ~ "Prioritized edges",
      TRUE ~ "Background"
    )
  )

p <- ggplot(
  plot_tbl,
  aes(x = recurrent_count, y = weight)
) +
  geom_point(
    data = plot_tbl %>% dplyr::filter(highlight_group == "Background"),
    color = "grey82",
    alpha = 0.7,
    size = 1.9
  ) +
  geom_point(
    data = plot_tbl %>% dplyr::filter(highlight_group == "Prioritized edges"),
    aes(color = direction),
    alpha = 0.95,
    size = 3.2
  ) +
  geom_point(
    data = plot_tbl %>% dplyr::filter(is_target),
    color = "black",
    fill = "#FFD54F",
    shape = 21,
    size = 4.2,
    stroke = 1.1
  ) +
  geom_text(
    data = label_tbl,
    aes(label = label),
    size = 3,
    nudge_y = 0.18,
    check_overlap = TRUE,
    color = "black",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(gain = "#B2182B", loss = "#2166AC"),
    name = "Direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Recurrent cancer count",
    y = "Recurrence × |Δ weight| (tumor vs normal)",
    title = "Prioritization of recurrent rewiring edges",
    subtitle = "Top annotated candidates are labeled; MMP9_SPP1 is highlighted as the exemplar edge"
  )

pdf(out_pdf, width = 9, height = 6.5)
print(p)
dev.off()

message("Saved: ", out_csv)
message("Saved: ", out_pdf)
