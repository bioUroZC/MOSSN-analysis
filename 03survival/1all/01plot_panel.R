rm(list = ls())

library(ggplot2)
library(patchwork)

base_dir <- "/proj/c.zihao/work1/03survival"
out_dir <- file.path(base_dir, "1all")

cancertypes <- c("ACC", "BLCA", "BRCA", "CHOL", "CRC", "GBM", "KIRC",
                 "LGG", "LIHC", "LUAD", "OV", "PAAD", "PRAD", "STAD")

method_levels <- c("MOSSN", "PPIXpress", "Patkar", "noCorr",
                   "noRWR", "RandomBackbone", "RawExpr", "NodeRWR")

custom_colors <- c(
  MOSSN          = "#b2182b",
  noRWR          = "#82a1c2",
  PPIXpress      = "#c9a1be",
  Patkar         = "#f8aa5d",
  RawExpr        = "#9ccdc9",
  NodeRWR        = "#87be81",
  noCorr         = "#ba9e8e",
  RandomBackbone = "#cfc8c5"
)

raw_data <- data.frame()
for (ct in cancertypes) {
  one <- read.csv(file.path(base_dir, ct, "03sur/ml_dataset.csv"), row.names = 1)
  one$CancerType <- ct
  raw_data <- rbind(raw_data, one)
}
raw_data$File <- factor(raw_data$File, levels = method_levels)

tumor_mean <- data.frame()
for (ct in cancertypes) {
  for (method in method_levels) {
    one <- raw_data[raw_data$CancerType == ct & raw_data$File == method, ]
    tumor_mean <- rbind(tumor_mean, data.frame(
      CancerType = ct,
      File = method,
      mean_C_index = mean(one$C_index),
      mean_tAUC = mean(one$Mean_tAUC)
    ))
  }
}
tumor_mean$CancerType <- factor(tumor_mean$CancerType, levels = cancertypes)
tumor_mean$File <- factor(tumor_mean$File, levels = method_levels)

out_tbl <- tumor_mean
out_tbl$mean_C_index <- round(out_tbl$mean_C_index, 5)
out_tbl$mean_tAUC <- round(out_tbl$mean_tAUC, 5)
write.csv(out_tbl, file.path(out_dir, "01tumor_mean.csv"), row.names = FALSE)

overall_means <- data.frame(
  File = method_levels,
  overall_mean_C_index = round(tapply(tumor_mean$mean_C_index, tumor_mean$File, mean), 5),
  overall_mean_tAUC = round(tapply(tumor_mean$mean_tAUC, tumor_mean$File, mean), 5)
)
print(overall_means[order(-overall_means$overall_mean_tAUC), ], row.names = FALSE)

small_theme <- theme_bw(base_size = 5) +
  theme(
    axis.text = element_text(size = 4),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )

make_bar <- function(y_var, y_label, title) {
  ggplot(tumor_mean, aes(x = CancerType, y = .data[[y_var]], fill = File)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = NULL, y = y_label, fill = "Method", title = title) +
    ylim(0, 0.9) +
    scale_fill_manual(values = custom_colors) +
    small_theme
}

p_a <- make_bar("mean_tAUC", "Mean timeAUC",
                "Mean timeAUC across 14 cancer types")
p_b <- make_bar("mean_C_index", "Mean C-index",
                "Mean C-index across 14 cancer types")

panel_ab <- (p_a / p_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 5, hjust = 0.5, margin = margin(b = 2)),
    plot.tag = element_text(size = 7, face = "bold")
  )

pdf(file.path(out_dir, "figure4ab.pdf"), height = 8 / 2.54, width = 8 / 2.54)
print(panel_ab)
dev.off()

rgb_val <- col2rgb(custom_colors) / 255
rgb_val <- rgb_val + (1 - rgb_val) * 0.3
point_colors <- setNames(rgb(rgb_val[1, ], rgb_val[2, ], rgb_val[3, ]), names(custom_colors))

make_box <- function(y_var, y_label, title) {
  ggplot(raw_data, aes(x = File, y = .data[[y_var]], fill = File)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.2) +
    geom_jitter(aes(color = File), width = 0.15, size = 0.3, alpha = 0.8) +
    labs(x = NULL, y = y_label, fill = "Method", title = title) +
    coord_cartesian(ylim = c(0.3, 1.0)) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = point_colors, guide = "none") +
    small_theme
}

p_c <- make_box("Mean_tAUC", "timeAUC",
                "Dataset-level distribution of timeAUC")
p_d <- make_box("C_index", "C-index",
                "Dataset-level distribution of C-index")

panel_cd <- (p_c | p_d) +
  plot_annotation(tag_levels = list(c("C", "D"))) &
  theme(
    legend.position = "none",
    plot.title = element_text(size = 5, hjust = 0.5, margin = margin(b = 2)),
    plot.title.position = "plot",
    plot.tag = element_text(size = 7, face = "bold")
  )

pdf(file.path(out_dir, "figure4cd.pdf"), height = 4 / 2.54, width = 8 / 2.54)
print(panel_cd)
dev.off()
