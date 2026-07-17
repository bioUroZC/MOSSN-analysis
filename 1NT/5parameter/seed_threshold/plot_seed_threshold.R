rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(ggplot2)
library(dplyr)

base_dir <- paste0(PROJ_ROOT, "/1NT/5parameter/seed_threshold")
auc_dir <- file.path(base_dir, "auc")
out_dir <- file.path(auc_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(auc_dir, "seed_threshold_metrics_summary.csv")
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)

quantile_levels <- c("q70", "q80", "q90", "q95")
quantile_labels <- c(q70 = "70th", q80 = "80th", q90 = "90th", q95 = "95th")

summary_df$seed_quantile <- factor(summary_df$seed_quantile, levels = quantile_levels)

metric_levels <- c("accuracy", "precision", "recall", "f1", "auc", "aupr")
metric_labels <- c(
  accuracy = "Accuracy",
  precision = "Precision",
  recall = "Recall",
  f1 = "F1",
  auc = "AUC",
  aupr = "AUPR"
)

metric_df <- bind_rows(lapply(metric_levels, function(metric_name) {
  data.frame(
    seed_quantile = summary_df$seed_quantile,
    metric = metric_name,
    value = summary_df[[metric_name]],
    stringsAsFactors = FALSE
  )
}))
metric_df$metric <- factor(metric_df$metric, levels = metric_levels, labels = metric_labels)

p1 <- ggplot(metric_df, aes(x = seed_quantile, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_x_discrete(labels = quantile_labels) +
  labs(x = "High-expression seed quantile", y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "seed_threshold_metrics.pdf"), p1, width = 8.4, height = 6.4)

metric_df_auc_accuracy <- metric_df %>%
  filter(metric %in% c("Accuracy", "AUC"))

p2 <- ggplot(metric_df_auc_accuracy, aes(x = seed_quantile, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_x_discrete(labels = quantile_labels) +
  labs(x = "High-expression seed quantile", y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "seed_threshold_accuracy_auc.pdf"), p2, width = 8.4, height = 3.4)
