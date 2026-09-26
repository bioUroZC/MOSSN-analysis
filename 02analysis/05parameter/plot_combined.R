rm(list = ls())

library(ggplot2)
library(patchwork)

base_dir <- "/proj/c.zihao/work1/02analysis/05parameter"

lambda_df <- read.csv(file.path(base_dir, "01lambda", "lambda_metrics_summary.csv"))
lambda_df$x <- as.numeric(sub("lam", "", lambda_df$lam)) / 10
lambda_df <- lambda_df[order(lambda_df$x), ]
lambda_df$x <- factor(format(lambda_df$x), levels = format(lambda_df$x))

restart_df <- read.csv(file.path(base_dir, "02restart", "restart_metrics_summary.csv"))
restart_df$x <- as.numeric(sub("alpha", "", restart_df$alpha)) / 10
restart_df <- restart_df[order(restart_df$x), ]
restart_df$x <- factor(format(restart_df$x), levels = format(restart_df$x))

seed_df <- read.csv(file.path(base_dir, "03seedThreshold", "seed_threshold_metrics_summary.csv"))
seed_df$x <- factor(
  seed_df$seed_quantile,
  levels = c("q70", "q80", "q90", "q95"),
  labels = c("70th", "80th", "90th", "95th")
)

plot_metric <- function(df, metric, x_lab, title) {
  df$y <- df[[metric]]

  ggplot(df, aes(x = x, y = y, group = 1)) +
    geom_line(linewidth = 0.4, color = "#3C5488") +
    geom_point(size = 0.8, color = "#3C5488") +
    labs(x = x_lab, y = "Metric value", title = title) +
    theme_bw(base_size = 6) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 5),
      axis.title = element_text(size = 5, color = "black"),
      axis.text = element_text(size = 5, color = "black"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

acc_title <- "Tumor-normal separation accuracy"
auc_title <- "Tumor-normal separation AUC"

p_a <- plot_metric(lambda_df, "accuracy", "Correction strength", acc_title)
p_b <- plot_metric(lambda_df, "auc", "Correction strength", auc_title)

restart_breaks <- scale_x_discrete(breaks = c("0.1", "0.3", "0.5", "0.7", "0.9"))
p_c <- plot_metric(restart_df, "accuracy", "Restart probability", acc_title) + restart_breaks
p_d <- plot_metric(restart_df, "auc", "Restart probability", auc_title) + restart_breaks
p_e <- plot_metric(seed_df, "accuracy", "High-expression seed quantile", acc_title)
p_f <- plot_metric(seed_df, "auc", "High-expression seed quantile", auc_title)

combined <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave(file.path(base_dir, "supfigure10.pdf"), combined,
       width = 10, height = 12, units = "cm")
