rm(list = ls())

library(ggplot2)
library(patchwork)

robust_dir <- "/proj/c.zihao/work1/02analysis/01robust"
out_dir <- file.path(robust_dir, "02merge")

retention_levels <- c(10, 30, 50, 70)
methods <- c("MOSSN", "Patkar", "PPIXpress", "Proteinarium", "SSN", "SWEET", "LIONESS")

method_colors <- c(
  MOSSN = "#b2182b", Patkar = "#f8aa5d", PPIXpress = "#c9a1be",
  SSN = "#82a1c2", SWEET = "#9ccdc9", Proteinarium = "#ba9e8e",
  LIONESS = "#87be81"
)

summary_df <- data.frame()
for (level in retention_levels) {
  one <- read.csv(file.path(robust_dir, level, "consistency", "consistency_df.csv"))
  summary_df <- rbind(summary_df, one)
}
summary_df <- summary_df[summary_df$method %in% methods, ]
summary_df$method <- factor(summary_df$method, levels = methods)
summary_df$retention_pct <- factor(summary_df$retention_pct, levels = retention_levels,
                                   labels = paste0(retention_levels, "%"))

write.csv(summary_df, file.path(out_dir, "robust_summary_all.csv"), row.names = FALSE)

make_barplot <- function(data, ncol, label_size = 3.2) {
  ggplot(data, aes(x = method, y = consistency, fill = method)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%.3f", consistency)), vjust = -0.35, size = label_size) +
    facet_wrap(~retention_pct, ncol = ncol) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(x = NULL, y = "Mean Spearman consistency") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )
}

small_theme <- theme(
  text = element_text(size = 5),
  axis.text = element_text(size = 4),
  strip.text = element_text(size = 5, margin = margin(1, 0, 1, 0)),
  strip.background = element_rect(fill = "grey95", linewidth = 0.3),
  panel.border = element_rect(fill = NA, linewidth = 0.3),
  panel.grid.major = element_line(linewidth = 0.2),
  axis.ticks = element_line(linewidth = 0.2),
  axis.ticks.length = unit(0.8, "mm"),
  plot.margin = margin(2, 2, 2, 2)
)

make_retention_panel <- function(pct, show_ylab) {
  title <- paste0("Spearman correlations between original and\n",
                  "subsampling-perturbed networks (", pct, " samples)")
  make_barplot(summary_df[summary_df$retention_pct == pct, ], ncol = 1,
               label_size = 1.1) +
    small_theme +
    labs(title = title, y = if (show_ylab) "Mean Spearman consistency" else NULL) +
    theme(
      strip.text = element_blank(),        
      strip.background = element_blank(),
      plot.title = element_text(size = 4, hjust = 0.5, margin = margin(b = 2)),
      plot.title.position = "plot"
    )
}

p_main <- make_retention_panel("10%", TRUE)
ggsave(file.path(out_dir, "figure3a.pdf"), p_main, width = 4, height = 4, units = "cm")

p_sup <- (make_retention_panel("30%", TRUE) | make_retention_panel("50%", FALSE) |
            make_retention_panel("70%", FALSE)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 7, face = "bold"))
ggsave(file.path(out_dir, "supfigure4.pdf"), p_sup, width = 12, height = 4, units = "cm")
