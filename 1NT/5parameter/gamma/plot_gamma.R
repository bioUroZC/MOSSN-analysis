rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(ggplot2)
library(dplyr)

base_dir <- paste0(PROJ_ROOT, "/1NT/5parameter/gamma")
auc_dir <- file.path(base_dir, "auc")
out_dir <- file.path(auc_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(auc_dir, "gamma_metrics_summary.csv")
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)

summary_df$gamma_value <- as.numeric(sub("g", "", summary_df$gamma)) / 10
gamma_breaks <- sort(unique(summary_df$gamma_value))
gamma_labels <- format(gamma_breaks, trim = TRUE, scientific = FALSE)

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
    gamma_value = summary_df$gamma_value,
    metric = metric_name,
    value = summary_df[[metric_name]],
    stringsAsFactors = FALSE
  )
}))
metric_df$gamma_factor <- factor(
  metric_df$gamma_value,
  levels = gamma_breaks,
  labels = gamma_labels
)
metric_df$metric <- factor(metric_df$metric, levels = metric_levels, labels = metric_labels)

p1 <- ggplot(metric_df, aes(x = gamma_factor, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = expression(gamma), y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "gamma_metrics.pdf"), p1, width = 8.4, height = 6.4)

metric_df_auc_accuracy <- metric_df %>%
  filter(metric %in% c("Accuracy", "AUC"))

p2 <- ggplot(metric_df_auc_accuracy, aes(x = gamma_factor, y = value, group = 1)) +
  geom_line(linewidth = 1, color = "#b22222") +
  geom_point(size = 2.2, color = "#b22222") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = expression(gamma), y = "Metric value") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "gamma_accuracy_auc.pdf"), p2, width = 8.4, height = 3.4)
