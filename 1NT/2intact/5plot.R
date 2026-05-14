rm(list = ls())

library(ggplot2)
library(dplyr)
library(cowplot)

base_dir <- "/proj/c.zihao/work1/1NT/2intact/1distance"
out_dir  <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_raw <- read.csv(file.path(base_dir, "results_Cluster.csv"),
                        stringsAsFactors = FALSE)

cancer_levels <- c("BLCA","BRCA","CRC","ESCA","HNSC","KIRC",
                   "LIHC","LUAD","LUSC","PRAD","STAD","Mean")

fraction_levels <- c("5%","10%","15%","20%")

method_levels <- c(
  "MOSSN_noPrior",
  "SSN","SWEET","Patkar","Proteinarium","PPIXpress","LIONESS"
)

method_colors <- c(
  MOSSN_noPrior   = "#67000d",
  SSN            = "#4c78a8",
  SWEET          = "#72b7b2",
  Patkar         = "#c9a227",
  Proteinarium   = "#9d755d",
  PPIXpress      = "#b279a2",
  LIONESS        = "#54a24b"
)

# ── data for line plots (all fractions, averaged across cancers) ───────────────
line_df <- cluster_raw %>%
  filter(method %in% method_levels) %>%
  mutate(
    method           = factor(method, levels = method_levels),
    feature_fraction = factor(feature_fraction, levels = fraction_levels)
  ) %>%
  group_by(method, feature_fraction) %>%
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_auc      = mean(auc,      na.rm = TRUE),
    .groups = "drop"
  )

# ── data for heatmaps (20% only) ──────────────────────────────────────────────
cluster_df <- cluster_raw %>%
  filter(method %in% method_levels, feature_fraction == "20%") %>%
  mutate(
    cancer = factor(cancer, levels = cancer_levels),
    method = factor(method, levels = method_levels)
  )

summary_df <- cluster_df %>%
  group_by(method, cancer) %>%
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_auc      = mean(auc,      na.rm = TRUE),
    .groups = "drop"
  )

mean_col <- summary_df %>%
  group_by(method) %>%
  summarise(
    mean_accuracy = mean(mean_accuracy, na.rm = TRUE),
    mean_auc      = mean(mean_auc,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cancer = "Mean")

summary_df <- bind_rows(summary_df, mean_col) %>%
  mutate(cancer = factor(cancer, levels = cancer_levels))

# ── shared theme ──────────────────────────────────────────────────────────────
heatmap_theme <- theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9.5,
                                    color = "grey20"),
    axis.text.y      = element_text(size = 9.5, color = "grey20"),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5,
                                    margin = margin(b = 6)),
    legend.position  = "right",
    legend.key.width  = unit(0.45, "cm"),
    legend.key.height = unit(1.6,  "cm"),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(8, 8, 8, 8)
  )

# diverging palette: blue (low) → white (0.5) → red (high)
div_colors <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
                "#f7f7f7",
                "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

make_heatmap <- function(df, value_col, title, legend_title) {
  n_cancer <- length(cancer_levels) - 1  # exclude "Mean"
  vals     <- df[[value_col]]
  lo       <- min(vals, na.rm = TRUE)
  hi       <- max(vals, na.rm = TRUE)
  df$txt_col <- ifelse(vals > 0.72, "white", "grey15")

  x_face <- c(rep("plain", n_cancer), "bold")

  ggplot(df, aes(x = cancer,
                 y = factor(method, levels = rev(method_levels)),
                 fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", .data[[value_col]]),
                  color = txt_col),
              size = 2.5, fontface = "plain") +
    # separator between last cancer and Mean column
    annotate("segment",
             x = n_cancer + 0.5, xend = n_cancer + 0.5,
             y = 0.5, yend = length(method_levels) + 0.5,
             color = "grey40", linewidth = 1.0) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors  = div_colors,
      values  = scales::rescale(c(0.5, 0.65, 0.8, 1.0), from = c(0.5, 1.0)),
      limits  = c(0.5, 1.0),
      oob     = scales::squish,
      name    = legend_title,
      guide   = guide_colorbar(frame.colour = "grey70", ticks.colour = "grey70")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = title) +
    heatmap_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5,
                                     color = "grey20", face = x_face))
}

p_acc <- make_heatmap(summary_df, "mean_accuracy", "Accuracy", "Accuracy")
p_auc <- make_heatmap(summary_df, "mean_auc",      "AUC",      "AUC")

ggsave(file.path(out_dir, "heatmap_accuracy.pdf"), p_acc,
       width = 8.5, height = 4.2)

# ── line plots ────────────────────────────────────────────────────────────────
line_theme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 9.5, color = "grey20"),
    axis.title         = element_text(size = 10),
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5,
                                      margin = margin(b = 6)),
    legend.position    = "right",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 8.5),
    legend.key.height  = unit(0.55, "cm"),
    legend.key.width   = unit(1.0,  "cm"),
    plot.margin        = margin(8, 8, 8, 8)
  )

make_line <- function(df, value_col, title, ylab) {
  ggplot(df, aes(x = feature_fraction, y = .data[[value_col]],
                 color = method, group = method)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2.0) +
    scale_color_manual(values = method_colors,
                       breaks = method_levels) +
    scale_y_continuous(limits = c(0.3, 1.01), breaks = seq(0.3, 1.0, 0.1)) +
    labs(x = "Feature fraction", y = ylab, title = title) +
    line_theme
}

p_line_auc <- make_line(line_df, "mean_auc",      "AUC",      "Mean AUC")
p_line_acc <- make_line(line_df, "mean_accuracy", "Accuracy", "Mean Accuracy")

# ── combined 3-panel figure: heatmap_auc (top) + line_acc + line_auc (bottom) ─
top_row <- plot_grid(p_line_acc, p_line_auc,
                     nrow = 1, labels = c("A", "B"),
                     label_size = 13)
combined <- plot_grid(top_row, p_auc,
                      nrow = 2, labels = c("", "C"),
                      label_size = 13,
                      rel_heights = c(1, 1))

ggsave(file.path(out_dir, "combined_auc_lines.pdf"), combined,
       width = 14, height = 8.4)

write.csv(summary_df, file.path(out_dir, "cluster_summary.csv"),
          row.names = FALSE)
