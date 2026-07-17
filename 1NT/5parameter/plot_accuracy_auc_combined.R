rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(ggplot2)
library(dplyr)
library(grid)

base_dir <- paste0(PROJ_ROOT, "/1NT/5parameter")
out_dir <- file.path(base_dir, "combined_auc", "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metric_levels <- c("accuracy", "auc")
metric_labels <- c(
  accuracy = "Accuracy",
  auc = "AUC"
)

build_metric_df <- function(summary_df, x_col, x_levels = NULL, x_labels = NULL) {
  metric_df <- bind_rows(lapply(metric_levels, function(metric_name) {
    data.frame(
      x_value = summary_df[[x_col]],
      metric = metric_name,
      value = summary_df[[metric_name]],
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(x_levels)) {
    metric_df$x_value <- factor(metric_df$x_value, levels = x_levels)
  }

  if (!is.null(x_labels)) {
    metric_df$x_label <- factor(
      metric_df$x_value,
      levels = if (is.null(x_levels)) unique(metric_df$x_value) else x_levels,
      labels = x_labels
    )
  } else {
    metric_df$x_label <- factor(metric_df$x_value)
  }

  metric_df$metric <- factor(metric_df$metric, levels = metric_levels, labels = metric_labels)
  metric_df
}

build_single_plot <- function(plot_df, metric_name, x_lab, panel_tag, parameter_name) {
  metric_display <- if (as.character(metric_name) == "AUC") {
    "LUAD paired tumor-normal PC1 ROC AUC"
  } else {
    "LUAD paired tumor-normal k-means accuracy"
  }

  ggplot(
    plot_df %>% filter(metric == metric_name),
    aes(x = x_label, y = value, group = 1)
  ) +
    geom_line(linewidth = 1, color = "#b22222") +
    geom_point(size = 2.2, color = "#b22222") +
    labs(
      x = x_lab,
      y = "Metric value",
      title = paste(
        panel_tag,
        paste("Effect of", parameter_name, "on", metric_display)
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold")
    )
}

gamma_file <- file.path(base_dir, "gamma", "auc", "gamma_metrics_summary.csv")
gamma_df <- read.csv(gamma_file, stringsAsFactors = FALSE)
gamma_df$gamma_value <- as.numeric(sub("g", "", gamma_df$gamma)) / 10
gamma_breaks <- sort(unique(gamma_df$gamma_value))
gamma_labels <- format(gamma_breaks, trim = TRUE, scientific = FALSE)
gamma_plot_df <- build_metric_df(
  summary_df = gamma_df,
  x_col = "gamma_value",
  x_levels = gamma_breaks,
  x_labels = gamma_labels
)

restart_file <- file.path(base_dir, "restart", "auc", "restart_metrics_summary.csv")
restart_df <- read.csv(restart_file, stringsAsFactors = FALSE)
restart_df$alpha_value <- as.numeric(sub("alpha", "", restart_df$alpha)) / 10
alpha_breaks <- sort(unique(restart_df$alpha_value))
alpha_labels <- format(alpha_breaks, trim = TRUE, scientific = FALSE)
restart_plot_df <- build_metric_df(
  summary_df = restart_df,
  x_col = "alpha_value",
  x_levels = alpha_breaks,
  x_labels = alpha_labels
)

seed_file <- file.path(base_dir, "seed_threshold", "auc", "seed_threshold_metrics_summary.csv")
seed_df <- read.csv(seed_file, stringsAsFactors = FALSE)
quantile_levels <- c("q70", "q80", "q90", "q95")
quantile_labels <- c("70th", "80th", "90th", "95th")
seed_plot_df <- build_metric_df(
  summary_df = seed_df,
  x_col = "seed_quantile",
  x_levels = quantile_levels,
  x_labels = quantile_labels
)

plot_list <- list(
  build_single_plot(gamma_plot_df, "Accuracy", expression(gamma), "A", "gamma"),
  build_single_plot(gamma_plot_df, "AUC", expression(gamma), "B", "gamma"),
  build_single_plot(restart_plot_df, "Accuracy", "Restart probability", "C", "restart probability"),
  build_single_plot(restart_plot_df, "AUC", "Restart probability", "D", "restart probability"),
  build_single_plot(seed_plot_df, "Accuracy", "High-expression seed quantile", "E", "seed threshold"),
  build_single_plot(seed_plot_df, "AUC", "High-expression seed quantile", "F", "seed threshold")
)

output_file <- file.path(out_dir, "combined_accuracy_auc.pdf")
pdf(output_file, width = 10, height = 12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))

for (i in seq_along(plot_list)) {
  row_idx <- ((i - 1) %/% 2) + 1
  col_idx <- ((i - 1) %% 2) + 1
  print(plot_list[[i]], vp = viewport(layout.pos.row = row_idx, layout.pos.col = col_idx))
}

dev.off()
