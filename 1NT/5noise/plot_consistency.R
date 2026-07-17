rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(ggplot2)
library(dplyr)

base_dir <- paste0(PROJ_ROOT, "/1NT/5noise")
levels <- c("05", "10", "15", "20")
level_labels <- c("05" = "5%", "10" = "10%", "15" = "15%", "20" = "20%")
out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c(
  "MOSSN_uniform", "Patkar", "PPIXpress", "Proteinarium",
  "SSN", "SWEET", "LIONESS"
)

method_colors <- c(
  MOSSN_uniform = "#67000d",
  SSN = "#4c78a8",
  SWEET = "#72b7b2",
  Patkar = "#c9a227",
  Proteinarium = "#9d755d",
  PPIXpress = "#b279a2",
  LIONESS = "#54a24b"
)

summary_list <- list()
sample_list <- list()

for (lvl in levels) {
  summary_file <- file.path(base_dir, lvl, "consistency", "consistency_df.csv")
  sample_file <- file.path(base_dir, lvl, "consistency", "sample_level_consistency.csv")

  if (file.exists(summary_file)) {
    df <- read.csv(summary_file, stringsAsFactors = FALSE)
    df$level <- lvl
    summary_list[[lvl]] <- df
  }

  if (file.exists(sample_file)) {
    df <- read.csv(sample_file, stringsAsFactors = FALSE)
    df$level <- lvl
    sample_list[[lvl]] <- df
  }
}

summary_df <- bind_rows(summary_list)
sample_df <- bind_rows(sample_list)

summary_df$method <- factor(summary_df$method, levels = method_levels)
summary_df$level <- factor(summary_df$level, levels = levels, labels = level_labels[levels])
sample_df$method <- factor(sample_df$method, levels = method_levels)
sample_df$level <- factor(sample_df$level, levels = levels, labels = level_labels[levels])

make_barplot <- function(data, ncol = 1) {
  ggplot(data, aes(x = method, y = consistency, fill = method)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%.3f", consistency)), vjust = -0.35, size = 3.2) +
    facet_wrap(~ level, ncol = ncol) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = "Mean Spearman consistency") +
    coord_cartesian(ylim = c(0, 1.02)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )
}

p_20 <- make_barplot(filter(summary_df, level == "20%"), ncol = 1)
ggsave(
  file.path(out_dir, "noise_barplot_20pct.pdf"),
  p_20,
  width = 6,
  height = 5
)

p_rest <- make_barplot(filter(summary_df, level != "20%"), ncol = 3)
ggsave(
  file.path(out_dir, "noise_barplot_5_10_15pct.pdf"),
  p_rest,
  width = 14,
  height = 5
)

write.csv(summary_df, file.path(out_dir, "noise_summary_all.csv"), row.names = FALSE)
write.csv(sample_df, file.path(out_dir, "noise_sample_level_all.csv"), row.names = FALSE)
