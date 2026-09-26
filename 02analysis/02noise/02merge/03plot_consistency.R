rm(list = ls())

library(ggplot2)
library(patchwork)

noise_dir <- "/proj/c.zihao/work1/02analysis/02noise"
out_dir <- file.path(noise_dir, "02merge")

noise_levels <- c(5, 10, 15, 20)
methods <- c("MOSSN", "Patkar", "PPIXpress", "Proteinarium", "SSN", "SWEET", "LIONESS")

method_colors <- c(
  MOSSN = "#b2182b", Patkar = "#f8aa5d", PPIXpress = "#c9a1be",
  SSN = "#82a1c2", SWEET = "#9ccdc9", Proteinarium = "#ba9e8e",
  LIONESS = "#87be81"
)

summary_df <- data.frame()
sample_df <- data.frame()
for (level in noise_levels) {
  level_dir <- file.path(noise_dir, sprintf("%02d", level), "consistency")
  summary_df <- rbind(summary_df, read.csv(file.path(level_dir, "consistency_df.csv")))
  sample_df <- rbind(sample_df, read.csv(file.path(level_dir, "sample_level_consistency.csv")))
}

summary_df <- summary_df[summary_df$method %in% methods, ]
summary_df$method <- factor(summary_df$method, levels = methods)
summary_df$noise_pct <- factor(summary_df$noise_pct, levels = noise_levels,
                               labels = paste0(noise_levels, "%"))
write.csv(summary_df, file.path(out_dir, "noise_summary_all.csv"), row.names = FALSE)

sample_df <- sample_df[sample_df$method %in% methods, ]
sample_df$method <- factor(sample_df$method, levels = methods)
sample_df$noise_pct <- factor(sample_df$noise_pct, levels = noise_levels,
                              labels = paste0(noise_levels, "%"))
write.csv(sample_df, file.path(out_dir, "noise_sample_level_all.csv"), row.names = FALSE)

make_barplot <- function(data, ncol, label_size = 3.2) {
  ggplot(data, aes(x = method, y = consistency, fill = method)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%.3f", consistency)), vjust = -0.35, size = label_size) +
    facet_wrap(~noise_pct, ncol = ncol) +
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

make_noise_panel <- function(pct, show_ylab) {
  title <- paste0("Spearman correlations between original and\n",
                  "noise-perturbed networks (", pct, " of the\n",
                  "gene-wise standard deviation)")
  make_barplot(summary_df[summary_df$noise_pct == pct, ], ncol = 1,
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

p_main <- make_noise_panel("20%", TRUE)
ggsave(file.path(out_dir, "figure3b.pdf"), p_main, width = 4, height = 4, units = "cm")

p_sup <- (make_noise_panel("5%", TRUE) | make_noise_panel("10%", FALSE) |
            make_noise_panel("15%", FALSE)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 7, face = "bold"))
ggsave(file.path(out_dir, "supfigure5.pdf"), p_sup, width = 12, height = 4, units = "cm")
