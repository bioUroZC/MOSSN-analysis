rm(list = ls())

library(ggplot2)
library(dplyr)
library(cowplot)

workflow_dir <- "/proj/c.zihao/work1/01NT/01string"
benchmark_dir <- file.path(workflow_dir, "02benchmark")
out_dir <- file.path(workflow_dir, "02benchmark", "03benchPlots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_raw <- read.csv(file.path(benchmark_dir, "02merge", "results_Cluster.csv"),
                        stringsAsFactors = FALSE)

cancer_levels <- c("BLCA","BRCA","CRC","ESCA","HNSC","KIRC",
                   "LIHC","LUAD","LUSC","PRAD","STAD","Mean")

fraction_levels <- c("5%","10%","15%","20%")

method_levels <- c(
  "MOSSN", "SSN","SWEET","LIONESS", "Patkar","Proteinarium","PPIXpress"
)

method_colors <- c(
  MOSSN          = "#b2182b",
  SSN            = "#82a1c2",
  SWEET          = "#9ccdc9",
  LIONESS        = "#87be81",
  Patkar         = "#f8aa5d",
  Proteinarium   = "#ba9e8e",
  PPIXpress      = "#c9a1be"

)

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

heatmap_theme <- theme_minimal(base_size = 5, base_family = "sans") +
  theme(
    axis.text.y      = element_text(size = 4, color = "grey20"),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(size = 6, hjust = 0.5,
                                    margin = margin(b = 2)),
    plot.title.position = "plot",  
    legend.position  = "right",
    legend.key.width  = unit(0.15, "cm"),
    legend.key.height = unit(0.5,  "cm"),  
    legend.title     = element_text(size = 4, face = "bold"),
    legend.text      = element_text(size = 3.5),
    legend.margin    = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    plot.margin      = margin(2, 2, 2, 2)
  )

div_colors <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
                "#f7f7f7",
                "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

make_heatmap <- function(df, value_col, title, legend_title) {
  n_cancer <- length(cancer_levels) - 1  
  vals     <- df[[value_col]]
  lo       <- min(vals, na.rm = TRUE)
  hi       <- max(vals, na.rm = TRUE)
  df$txt_col <- ifelse(vals >= 0.90 | vals <= 0.55, "white", "grey15")

  x_face <- c(rep("plain", n_cancer), "bold")

  ggplot(df, aes(x = cancer,
                 y = factor(method, levels = rev(method_levels)),
                 fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3f", .data[[value_col]]),
                  color = txt_col),
              size = 1.3, fontface = "plain") +
    
    annotate("segment",
             x = n_cancer + 0.5, xend = n_cancer + 0.5,
             y = 0.5, yend = length(method_levels) + 0.5,
             color = "grey40", linewidth = 0.4) +
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4,
                                     color = "grey20", face = x_face),
          legend.position = "none")
}

p_acc <- make_heatmap(summary_df, "mean_accuracy",
                      "Tumor-normal separation accuracy in 20% top-ranked interactions",
                      "Accuracy")
p_auc <- make_heatmap(summary_df, "mean_auc",
                      "Tumor-normal separation AUC in 20% top-ranked interactions",
                      "AUC")

ggsave(file.path(out_dir, "figure2c.pdf"), p_acc,
       width = 8, height = 4, units = "cm")
ggsave(file.path(out_dir, "figure2d.pdf"), p_auc,
       width = 8, height = 4, units = "cm")

line_theme <- theme_bw(base_size = 6, base_family = "sans") +
  theme(
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 5, color = "grey20"),
    axis.title         = element_text(size = 6),
    plot.title         = element_text(size = 6, hjust = 0.5,
                                      margin = margin(b = 2)),
    legend.position    = "right",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 5),
    legend.key.height  = unit(0.3, "cm"),
    legend.key.width   = unit(0.3,  "cm"),
    legend.margin      = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    plot.margin        = margin(2, 2, 2, 2)
  )

make_line <- function(df, value_col, title, ylab) {
  ggplot(df, aes(x = feature_fraction, y = .data[[value_col]],
                 color = method, group = method)) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 0.7) +
    scale_color_manual(values = method_colors,
                       breaks = method_levels) +
    scale_y_continuous(limits = c(0.3, 1.01), breaks = seq(0.3, 1.0, 0.1)) +
    labs(x = "Feature fraction", y = ylab, title = title) +
    line_theme
}

p_line_auc <- make_line(line_df, "mean_auc",
                        "Tumor-normal separation AUC", "Mean AUC")
p_line_acc <- make_line(line_df, "mean_accuracy",
                        "Tumor-normal separation accuracy", "Mean Accuracy")

ggsave(file.path(out_dir, "figure2b.pdf"), p_line_auc,
       width = 6, height = 4, units = "cm")
ggsave(file.path(out_dir, "figure2A.pdf"), p_line_acc,
       width = 6, height = 4, units = "cm")

write.csv(summary_df, file.path(out_dir, "cluster_summary.csv"),
          row.names = FALSE)
