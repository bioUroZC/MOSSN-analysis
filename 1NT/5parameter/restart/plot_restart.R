rm(list = ls())

library(ggplot2)
library(dplyr)

base_dir <- "/proj/c.zihao/work1/1NT/5parameter/restart"
auc_dir <- file.path(base_dir, "auc")
out_dir <- file.path(auc_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(auc_dir, "restart_metrics_summary.csv")
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)

summary_df$alpha_value <- as.numeric(sub("alpha", "", summary_df$alpha)) / 10
alpha_breaks <- sort(unique(summary_df$alpha_value))
alpha_labels <- format(alpha_breaks, trim = TRUE, scientific = FALSE)

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
    alpha_value = summary_df$alpha_value,
    metric = metric_name,
    value = summary_df[[metric_name]],
    stringsAsFactors = FALSE
  )
}))
metric_df$alpha_factor <- factor(
  metric_df$alpha_value,
  levels = alpha_breaks,
  labels = alpha_labels
)
metric_df$metric <- factor(metric_df$metric, levels = metric_levels, labels = metric_labels)

p1 <- ggplot(metric_df, aes(x = alpha_factor, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = "Restart probability", y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "restart_metrics.pdf"), p1, width = 8.4, height = 6.4)

metric_df_auc_accuracy <- metric_df %>%
  filter(metric %in% c("Accuracy", "AUC"))

p2 <- ggplot(metric_df_auc_accuracy, aes(x = alpha_factor, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = "Restart probability", y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "restart_accuracy_auc.pdf"), p2, width = 8.4, height = 3.4)
