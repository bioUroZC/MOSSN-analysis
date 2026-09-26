rm(list = ls())

library(ggplot2)
library(patchwork)
library(scales)

base_dir <- "/proj/c.zihao/work1/05atlas"
module_dir <- file.path(base_dir, "module_results")

cancers <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
             "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

gained_color <- "#d6604d"
lost_color <- "#4393c3"
base_family <- "Helvetica"

theme_atlas <- theme_bw(base_size = 5, base_family = base_family) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "#222222", size = 3),
    axis.title = element_text(size = 4.3),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.6, "mm"),
    plot.title = element_text(size = 3),
    plot.title.position = "plot",
    plot.subtitle = element_text(size = 4.2, color = "#444444"),
    plot.margin = margin(2, 2, 2, 2, "pt"),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key.size = unit(1.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.spacing = unit(1, "mm"),
    panel.border = element_rect(color = "#555555", linewidth = 0.25)
  )

atlas_all <- read.csv(file.path(base_dir, "atlas_all.csv"), row.names = 1)
atlas_all <- atlas_all[atlas_all$cancer %in% cancers, ]
recurrent_tbl <- read.csv(file.path(base_dir, "universal_recurrent_links.csv"))
module_summary <- read.csv(file.path(module_dir, "module_summary.csv"))
module_annotations <- read.csv(file.path(module_dir, "module_annotation_table.csv"))

cancer_summary <- data.frame()
for (cancer_name in cancers) {
  one <- atlas_all[atlas_all$cancer == cancer_name, ]
  cancer_summary <- rbind(cancer_summary, data.frame(
    cancer = cancer_name,
    n_pairs = one$n_pairs[1],
    tau_used = one$tau_used[1],
    n_tested = nrow(one),
    n_gained = sum(one$direction == "gained"),
    n_lost = sum(one$direction == "lost")
  ))
}
write.csv(cancer_summary, file.path(base_dir, "atlas_cancer_summary.csv"), row.names = FALSE)

plot_a_df <- data.frame(
  cancer = factor(rep(cancers, each = 2), levels = cancers),
  direction = rep(c("gained", "lost"), times = length(cancers)),
  n_links = as.vector(rbind(cancer_summary$n_gained, cancer_summary$n_lost))
)
plot_a_df$y <- ifelse(plot_a_df$direction == "lost", -plot_a_df$n_links, plot_a_df$n_links)

offset_a <- max(plot_a_df$n_links) * 0.025
p_a <- ggplot(plot_a_df, aes(cancer, y, fill = direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_text(
    aes(label = n_links, y = ifelse(direction == "gained", y + offset_a, y - offset_a),
        vjust = ifelse(direction == "gained", 0, 1)),
    size = 1.1, family = base_family
  ) +
  scale_fill_manual(values = c(gained = gained_color, lost = lost_color)) +
  labs(title = "Pan-cancer gained/lost interactions", x = NULL, y = "Number of links") +
  theme_atlas +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

jaccard_tbl <- data.frame()
for (direction in c("gained", "lost")) {
  for (i in seq_along(cancers)) {
    for (j in seq_along(cancers)) {
      set_i <- atlas_all$link[atlas_all$cancer == cancers[i] & atlas_all$direction == direction]
      set_j <- atlas_all$link[atlas_all$cancer == cancers[j] & atlas_all$direction == direction]
      union_n <- length(union(set_i, set_j))
      intersection_n <- length(intersect(set_i, set_j))
      jaccard_tbl <- rbind(jaccard_tbl, data.frame(
        cancer_x = cancers[j],
        cancer_y = cancers[i],
        x_idx = j,
        y_idx = i,
        direction = direction,
        intersection_n = intersection_n,
        union_n = union_n,
        jaccard = ifelse(union_n > 0, intersection_n / union_n, NA)
      ))
    }
  }
}
write.csv(jaccard_tbl, file.path(base_dir, "atlas_pairwise_jaccard.csv"), row.names = FALSE)

plot_b_df <- expand.grid(cancer_x = cancers, cancer_y = cancers, stringsAsFactors = FALSE)
plot_b_df$x_idx <- match(plot_b_df$cancer_x, cancers)
plot_b_df$y_idx <- match(plot_b_df$cancer_y, cancers)

key <- paste(plot_b_df$cancer_x, plot_b_df$cancer_y)
gained_tbl <- jaccard_tbl[jaccard_tbl$direction == "gained", ]
lost_tbl <- jaccard_tbl[jaccard_tbl$direction == "lost", ]
gained_jaccard <- gained_tbl$jaccard[match(key, paste(gained_tbl$cancer_x, gained_tbl$cancer_y))]
lost_jaccard <- lost_tbl$jaccard[match(key, paste(lost_tbl$cancer_x, lost_tbl$cancer_y))]

plot_b_df$fill_value <- NA
upper <- plot_b_df$x_idx > plot_b_df$y_idx
lower <- plot_b_df$x_idx < plot_b_df$y_idx
plot_b_df$fill_value[upper] <- gained_jaccard[upper]
plot_b_df$fill_value[lower] <- -lost_jaccard[lower]
plot_b_df$label <- percent(abs(plot_b_df$fill_value), accuracy = 1)
plot_b_df$cancer_x <- factor(plot_b_df$cancer_x, levels = cancers)
plot_b_df$cancer_y <- factor(plot_b_df$cancer_y, levels = rev(cancers))
plot_b_df <- plot_b_df[upper | lower, ]

max_jaccard <- max(abs(plot_b_df$fill_value))
p_b <- ggplot(plot_b_df, aes(cancer_x, cancer_y)) +
  geom_tile(aes(fill = fill_value), color = "white", linewidth = 0.15) +
  geom_text(aes(label = label), size = 0.85, family = base_family) +
  scale_fill_gradient2(
    low = lost_color, mid = "white", high = gained_color, midpoint = 0,
    limits = c(-max_jaccard, max_jaccard),
    breaks = c(-max_jaccard, 0, max_jaccard),
    labels = c("Lost", "", "Gained"),
    name = "Jaccard"
  ) +
  coord_fixed() +
  labs(title = "Recurrent gained and lost interactions", x = NULL, y = NULL) +
  theme_atlas +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 3),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

recurrent_class <- ifelse(recurrent_tbl$recurrent_class == "recurrently_gained", "Gained", "Lost")
plot_c_df <- as.data.frame(table(
  recurrent_count = factor(recurrent_tbl$recurrent_count, levels = 7:11),
  class = recurrent_class
))
names(plot_c_df)[3] <- "n_links"
plot_c_df <- plot_c_df[plot_c_df$n_links > 0, ]
plot_c_df <- plot_c_df[order(plot_c_df$recurrent_count, plot_c_df$class), ]

p_c <- ggplot(plot_c_df, aes(recurrent_count, n_links, fill = class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n_links), position = position_stack(vjust = 0.5),
            size = 1.3, family = base_family) +
  scale_fill_manual(values = c(Gained = gained_color, Lost = lost_color)) +
  labs(
    title = "Cross-cancer recurrence of interaction changes",
    x = "Number of cancer types",
    y = "Number of recurrent interactions",
    fill = "Class"
  ) +
  theme_atlas +
  theme(legend.position = "top")

heatmap_edges <- recurrent_tbl
heatmap_edges$abs_median_delta <- abs(heatmap_edges$median_delta_across_cancers)
heatmap_edges <- heatmap_edges[order(heatmap_edges$recurrent_class,
                                     -heatmap_edges$recurrent_count,
                                     -heatmap_edges$abs_median_delta,
                                     heatmap_edges$link, method = "radix"), ]
heatmap_edges <- rbind(
  head(heatmap_edges[heatmap_edges$recurrent_class == "recurrently_gained", ], 12),
  head(heatmap_edges[heatmap_edges$recurrent_class == "recurrently_lost", ], 12)
)
heatmap_edges$class <- ifelse(heatmap_edges$recurrent_class == "recurrently_gained",
                              "Recurrently gained", "Recurrently lost")
heatmap_edges$edge_order <- seq_len(nrow(heatmap_edges))
write.csv(heatmap_edges, file.path(base_dir, "atlas_heatmap_edges.csv"), row.names = FALSE)

plot_d_df <- atlas_all[atlas_all$link %in% heatmap_edges$link,
                       c("link", "cancer", "direction", "delta_median")]
plot_d_df <- plot_d_df[order(match(plot_d_df$link, heatmap_edges$link)), ]
plot_d_df$class <- heatmap_edges$class[match(plot_d_df$link, heatmap_edges$link)]
plot_d_df$plot_delta <- ifelse(plot_d_df$direction %in% c("gained", "lost"), plot_d_df$delta_median, NA)
plot_d_df$link <- factor(plot_d_df$link, levels = rev(heatmap_edges$link))
plot_d_df$cancer <- factor(plot_d_df$cancer, levels = cancers)

heat_limit <- quantile(abs(plot_d_df$plot_delta), 0.98, na.rm = TRUE)
p_d <- ggplot(plot_d_df, aes(cancer, link, fill = plot_delta)) +
  geom_tile(color = "white", linewidth = 0.12) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = lost_color, mid = "white", high = gained_color, midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = squish,
    na.value = "#eeeeee", name = "Median difference"
  ) +
  labs(title = "Most recurrent interaction perturbations", x = NULL, y = NULL) +
  theme_atlas +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
    axis.text.y = element_text(size = 3),
    strip.text = element_text(size = 3.2, margin = margin(1, 1, 1, 1, "pt")),
    legend.position = "right"
  )

label_e <- data.frame()
for (direction in c("gained", "lost")) {
  one <- module_summary[module_summary$direction == direction, ]
  one <- one[order(-one$total_weight), ]
  label_e <- rbind(label_e, head(one, 3))
}

p_e <- ggplot(module_summary, aes(n_genes, total_weight)) +
  geom_point(aes(size = n_edges, color = direction), alpha = 0.78, stroke = 0.2) +
  geom_text(data = label_e, aes(label = module_id),
            size = 1.2, nudge_y = 0.12, check_overlap = TRUE,
            family = base_family, show.legend = FALSE) +
  scale_y_log10(labels = label_number()) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  scale_color_manual(values = c(gained = gained_color, lost = lost_color),
                     labels = c(gained = "Gained", lost = "Lost")) +
  scale_size_continuous(range = c(0.5, 2.8)) +
  labs(
    title = "Organization of recurrent rewiring into modules",
    x = "Number of genes",
    y = "Total module weight",
    color = "Direction",
    size = "Edges"
  ) +
  theme_atlas +
  theme(legend.position = "right")

plot_f_df <- module_annotations[module_annotations$annotation_status == "ok", ]
plot_f_df$best_fdr <- ifelse(plot_f_df$best_source_hk == "Hallmark",
                             plot_f_df$hallmark_fdr, plot_f_df$kegg_fdr)
plot_f_df$annotation_score <- -log10(pmax(plot_f_df$best_fdr, 1e-300))

top_f <- data.frame()
for (direction in c("gained", "lost")) {
  one <- plot_f_df[plot_f_df$direction == direction, ]
  one <- one[order(-one$annotation_score), ]
  top_f <- rbind(top_f, head(one, 5))
}
plot_f_df <- top_f[order(top_f$direction, top_f$annotation_score, method = "radix"), ]
plot_f_df$annotation_label <- paste0(plot_f_df$module_id, " | ", plot_f_df$concise_label)
plot_f_df$annotation_label <- factor(plot_f_df$annotation_label,
                                     levels = unique(plot_f_df$annotation_label))

p_f <- ggplot(plot_f_df, aes(annotation_score, annotation_label, fill = direction)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(gained = gained_color, lost = lost_color),
                    labels = c(gained = "Gained", lost = "Lost")) +
  labs(
    title = "Functional annotation of recurrent modules",
    x = "-log10(FDR)",
    y = NULL,
    fill = "Direction"
  ) +
  theme_atlas +
  theme(
    axis.text.y = element_text(size = 3),
    legend.position = "top",
    legend.justification = "right"
  )

metrics_tbl <- data.frame(
  metric = c("cancer_types", "tested_interactions", "recurrent_interactions",
             "recurrently_gained", "recurrently_lost", "gained_modules",
             "lost_modules", "annotated_modules"),
  value = c(
    length(cancers),
    length(unique(atlas_all$link)),
    nrow(recurrent_tbl),
    sum(recurrent_tbl$recurrent_class == "recurrently_gained"),
    sum(recurrent_tbl$recurrent_class == "recurrently_lost"),
    sum(module_summary$direction == "gained"),
    sum(module_summary$direction == "lost"),
    sum(module_annotations$annotation_status == "ok")
  )
)
write.csv(metrics_tbl, file.path(base_dir, "atlas_summary_metrics.csv"), row.names = FALSE)

top_row <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1, 0.9))
bottom_row <- (p_d | p_e | p_f) + plot_layout(widths = c(1, 0.9, 1.1))
p_final <- (top_row / bottom_row) +
  plot_layout(heights = c(0.9, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 7, family = base_family))

ggsave(file.path(base_dir, "Figure6.pdf"), p_final, width = 14, height = 8, units = "cm")
