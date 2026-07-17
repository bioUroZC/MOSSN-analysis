rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- paste0(PROJ_ROOT, "/1NT/9analysis")
module_dir <- file.path(base_dir, "module_results")
setwd(base_dir)

out_pdf <- file.path(base_dir, "figure6_perturbation_atlas.pdf")

cancers <- c(
  "BLCA", "BRCA", "CRC", "ESCA", "HNSC",
  "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

gain_color <- "#b2182b"
loss_color <- "#2166ac"
base_family <- "Helvetica"

theme_atlas <- theme_bw(base_size = 10.5, base_family = base_family) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "#222222"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 8.5, color = "#444444"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "#555555", linewidth = 0.5)
  )

atlas_all <- read.csv(
  file.path(base_dir, "atlas_all.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (names(atlas_all)[1] %in% c("", "X")) {
  atlas_all <- atlas_all[, -1, drop = FALSE]
}

recurrent_tbl <- read.csv(
  file.path(base_dir, "universal_recurrent_links.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
module_summary <- read.csv(
  file.path(module_dir, "module_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
module_annotations <- read.csv(
  file.path(module_dir, "module_annotation_table.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# A. Cancer-specific perturbation burden and a compact audit table.
cancer_summary <- atlas_all %>%
  filter(cancer %in% cancers) %>%
  group_by(cancer) %>%
  summarise(
    n_pairs = first(n_pairs),
    tau_used = first(tau_used),
    n_tested = n(),
    n_gain = sum(direction == "gain", na.rm = TRUE),
    n_loss = sum(direction == "loss", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cancer = factor(cancer, levels = cancers)) %>%
  arrange(cancer)

write.csv(cancer_summary, "atlas_cancer_summary.csv", row.names = FALSE)

plot_a_df <- cancer_summary %>%
  select(cancer, n_gain, n_loss) %>%
  pivot_longer(c(n_gain, n_loss), names_to = "direction", values_to = "n_links") %>%
  mutate(
    direction = recode(direction, n_gain = "gain", n_loss = "loss"),
    y = if_else(direction == "loss", -n_links, n_links)
  )

offset_a <- max(plot_a_df$n_links) * 0.025
p_a <- ggplot(plot_a_df, aes(cancer, y, fill = direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_text(
    aes(label = n_links, y = ifelse(direction == "gain", y + offset_a, y - offset_a)),
    size = 2.55, family = base_family
  ) +
  scale_fill_manual(values = c(gain = gain_color, loss = loss_color)) +
  labs(
    title = "Tumor-normal interaction perturbation burden",
    x = NULL,
    y = "Number of interactions\n(loss shown as negative)"
  ) +
  theme_atlas +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# B. Pairwise overlap normalized by the union size (Jaccard index).
sig_sets <- atlas_all %>%
  filter(cancer %in% cancers, direction %in% c("gain", "loss")) %>%
  select(cancer, direction, link)

jaccard_rows <- list()
for (direction_label in c("gain", "loss")) {
  for (i in seq_along(cancers)) {
    for (j in seq_along(cancers)) {
      ca_i <- cancers[i]
      ca_j <- cancers[j]
      set_i <- sig_sets$link[sig_sets$cancer == ca_i & sig_sets$direction == direction_label]
      set_j <- sig_sets$link[sig_sets$cancer == ca_j & sig_sets$direction == direction_label]
      union_n <- length(union(set_i, set_j))
      intersection_n <- length(intersect(set_i, set_j))
      jaccard_rows[[length(jaccard_rows) + 1]] <- data.frame(
        cancer_x = ca_j,
        cancer_y = ca_i,
        x_idx = j,
        y_idx = i,
        direction = direction_label,
        intersection_n = intersection_n,
        union_n = union_n,
        jaccard = ifelse(union_n > 0, intersection_n / union_n, NA_real_)
      )
    }
  }
}
jaccard_tbl <- bind_rows(jaccard_rows)
write.csv(jaccard_tbl, "atlas_pairwise_jaccard.csv", row.names = FALSE)

plot_b_df <- expand.grid(
  cancer_x = cancers,
  cancer_y = cancers,
  stringsAsFactors = FALSE
) %>%
  mutate(x_idx = match(cancer_x, cancers), y_idx = match(cancer_y, cancers)) %>%
  left_join(
    jaccard_tbl %>% filter(direction == "gain") %>%
      select(cancer_x, cancer_y, gain_jaccard = jaccard),
    by = c("cancer_x", "cancer_y")
  ) %>%
  left_join(
    jaccard_tbl %>% filter(direction == "loss") %>%
      select(cancer_x, cancer_y, loss_jaccard = jaccard),
    by = c("cancer_x", "cancer_y")
  ) %>%
  mutate(
    triangle = case_when(x_idx > y_idx ~ "gain", x_idx < y_idx ~ "loss", TRUE ~ "diag"),
    fill_value = case_when(
      triangle == "gain" ~ gain_jaccard,
      triangle == "loss" ~ -loss_jaccard,
      TRUE ~ NA_real_
    ),
    label = ifelse(triangle == "diag", NA_character_, percent(abs(fill_value), accuracy = 1)),
    cancer_x = factor(cancer_x, levels = cancers),
    cancer_y = factor(cancer_y, levels = rev(cancers))
  )

max_jaccard <- max(abs(plot_b_df$fill_value), na.rm = TRUE)
p_b <- ggplot(plot_b_df %>% filter(triangle != "diag"), aes(cancer_x, cancer_y)) +
  geom_tile(aes(fill = fill_value), color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 2.15, family = base_family) +
  scale_fill_gradient2(
    low = loss_color, mid = "white", high = gain_color, midpoint = 0,
    limits = c(-max_jaccard, max_jaccard),
    breaks = c(-max_jaccard, 0, max_jaccard),
    labels = c(
      paste0("Lost ", percent(max_jaccard, accuracy = 1)),
      "0",
      paste0("Gained ", percent(max_jaccard, accuracy = 1))
    ),
    name = "Jaccard"
  ) +
  coord_fixed() +
  labs(
    title = "Normalized overlap between cancer types",
    subtitle = "Upper: gained interactions; lower: lost interactions",
    x = NULL, y = NULL
  ) +
  theme_atlas +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

# C. Recurrence distribution supporting the seven-cancer definition.
plot_c_df <- recurrent_tbl %>%
  mutate(
    class = recode(
      recurrent_class,
      recurrently_gained = "Gained",
      recurrently_lost = "Lost"
    ),
    recurrent_count = factor(recurrent_count, levels = 7:11)
  ) %>%
  count(recurrent_count, class, name = "n_links")

p_c <- ggplot(plot_c_df, aes(recurrent_count, n_links, fill = class)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = n_links), position = position_stack(vjust = 0.5),
    size = 2.7, family = base_family
  ) +
  scale_fill_manual(values = c(Gained = gain_color, Lost = loss_color)) +
  labs(
    title = "Cross-cancer recurrence of interaction changes",
    x = "Number of cancer types",
    y = "Number of recurrent interactions",
    fill = "Class"
  ) +
  theme_atlas +
  theme(legend.position = "top")

# D. An objective heatmap: top 12 edges per direction by recurrence, then effect size.
heatmap_edges <- recurrent_tbl %>%
  mutate(abs_median_delta = abs(median_delta_across_cancers)) %>%
  group_by(recurrent_class) %>%
  arrange(desc(recurrent_count), desc(abs_median_delta), link, .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(
    class = recode(
      recurrent_class,
      recurrently_gained = "Recurrently gained",
      recurrently_lost = "Recurrently lost"
    ),
    edge_order = row_number()
  )
write.csv(heatmap_edges, "atlas_heatmap_edges.csv", row.names = FALSE)

plot_d_df <- heatmap_edges %>%
  select(link, class, edge_order) %>%
  left_join(
    atlas_all %>% select(link, cancer, direction, delta_median),
    by = "link"
  ) %>%
  filter(cancer %in% cancers) %>%
  mutate(
    plot_delta = if_else(direction %in% c("gain", "loss"), delta_median, NA_real_),
    link = factor(link, levels = rev(unique(heatmap_edges$link))),
    cancer = factor(cancer, levels = cancers)
  )

heat_limit <- quantile(abs(plot_d_df$plot_delta), 0.98, na.rm = TRUE)
p_d <- ggplot(plot_d_df, aes(cancer, link, fill = plot_delta)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = loss_color, mid = "white", high = gain_color, midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = squish,
    na.value = "#eeeeee", name = "Median difference"
  ) +
  labs(
    title = "Most recurrent interaction perturbations",
    subtitle = "Top 12 per direction by recurrence and |median difference|; grey indicates no significant change",
    x = NULL, y = NULL
  ) +
  theme_atlas +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6.5),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "right"
  )

# E. Module organization without selecting an exemplar edge.
label_e <- module_summary %>%
  group_by(direction) %>%
  slice_max(total_weight, n = 3, with_ties = FALSE) %>%
  ungroup()

p_e <- ggplot(module_summary, aes(n_genes, total_weight)) +
  geom_point(aes(size = n_edges, color = direction), alpha = 0.78) +
  geom_text(
    data = label_e, aes(label = module_id),
    size = 2.5, nudge_y = 0.09, check_overlap = TRUE,
    family = base_family, show.legend = FALSE
  ) +
  scale_y_log10(labels = label_number()) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  scale_color_manual(values = c(gain = gain_color, loss = loss_color)) +
  scale_size_continuous(range = c(1.8, 8)) +
  labs(
    title = "Organization of recurrent rewiring into modules",
    x = "Number of genes",
    y = "Total module weight (log scale)",
    color = "Direction",
    size = "Edges"
  ) +
  theme_atlas +
  theme(legend.position = "right")

# F. Strongest functional annotations per direction.
plot_f_df <- module_annotations %>%
  filter(annotation_status == "ok") %>%
  mutate(
    best_fdr = case_when(
      best_source_hk == "Hallmark" ~ hallmark_fdr,
      best_source_hk == "KEGG" ~ kegg_fdr,
      TRUE ~ NA_real_
    ),
    annotation_score = -log10(pmax(best_fdr, 1e-300))
  ) %>%
  group_by(direction) %>%
  slice_max(annotation_score, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(direction, annotation_score) %>%
  mutate(
    annotation_label = paste0(module_id, " | ", concise_label),
    annotation_label = factor(annotation_label, levels = unique(annotation_label))
  )

p_f <- ggplot(plot_f_df, aes(annotation_score, annotation_label, fill = direction)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(gain = gain_color, loss = loss_color)) +
  labs(
    title = "Functional annotation of recurrent modules",
    x = expression(-log[10](FDR)),
    y = NULL,
    fill = "Direction"
  ) +
  theme_atlas +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "top"
  )

metrics_tbl <- data.frame(
  metric = c(
    "cancer_types", "tested_interactions", "recurrent_interactions",
    "recurrently_gained", "recurrently_lost", "gain_modules",
    "loss_modules", "annotated_modules"
  ),
  value = c(
    length(cancers), length(unique(atlas_all$link)), nrow(recurrent_tbl),
    sum(recurrent_tbl$recurrent_class == "recurrently_gained"),
    sum(recurrent_tbl$recurrent_class == "recurrently_lost"),
    sum(module_summary$direction == "gain"),
    sum(module_summary$direction == "loss"),
    sum(module_annotations$annotation_status == "ok")
  )
)
write.csv(metrics_tbl, "atlas_summary_metrics.csv", row.names = FALSE)

top_row <- p_a | p_b | p_c
bottom_row <- p_d | p_e | p_f
p_final <- (top_row / bottom_row) +
  plot_layout(heights = c(0.9, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14, family = base_family))

ggsave(out_pdf, p_final, width = 18, height = 10.5)

message("Saved: ", out_pdf)
message("Saved atlas audit tables in: ", base_dir)
