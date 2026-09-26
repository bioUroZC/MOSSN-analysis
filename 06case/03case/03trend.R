rm(list = ls())

library(ggplot2)

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")

module_order <- c("IFNG_response", "Antigen_presentation",
                  "Cytotoxicity", "Checkpoint_neighborhood")
module_label <- c(
  IFNG_response = "IFN-gamma response",
  Antigen_presentation = "Antigen presentation",
  Cytotoxicity = "Cytotoxicity",
  Checkpoint_neighborhood = "Checkpoint neighborhood"
)
response_levels <- c("CR", "PR", "SD", "PD")

score_long <- read.csv(file.path(out_dir, "module_scores_long.csv"))
score_wide <- read.csv(file.path(out_dir, "module_scores_wide.csv"))

idx <- match(score_long$Sample, score_wide$Sample)
score_long$Response <- score_wide$Response[idx]
score_long$response_ordinal <- score_wide$response_ordinal[idx]
score_long$Response <- factor(score_long$Response, levels = response_levels)

trend_df <- data.frame()
for (module_name in module_order) {
  sub <- score_long[score_long$module == module_name, ]

  for (score_type in c("Expression", "Interaction")) {
    if (score_type == "Expression") {
      score <- sub$expression_score
    } else {
      score <- sub$interaction_score
    }

    kw <- kruskal.test(score ~ sub$Response)
    sp <- cor.test(score, sub$response_ordinal, method = "spearman", exact = FALSE)

    trend_df <- rbind(trend_df, data.frame(
      module = module_name,
      score_type = score_type,
      kw_p_value = kw$p.value,
      spearman_rho = as.numeric(sp$estimate),
      spearman_p_value = sp$p.value
    ))
  }
}

trend_df$kw_fdr <- p.adjust(trend_df$kw_p_value, method = "BH")
trend_df$spearman_fdr <- p.adjust(trend_df$spearman_p_value, method = "BH")

plot_df <- rbind(
  data.frame(module = score_long$module, Response = score_long$Response,
             score_type = "Expression", score = score_long$expression_score),
  data.frame(module = score_long$module, Response = score_long$Response,
             score_type = "Interaction", score = score_long$interaction_score)
)
plot_df$module_label <- factor(module_label[plot_df$module], levels = module_label[module_order])
plot_df$score_type <- factor(plot_df$score_type, levels = c("Expression", "Interaction"))

label_df <- data.frame(
  module_label = factor(module_label[trend_df$module], levels = module_label[module_order]),
  score_type = factor(trend_df$score_type, levels = c("Expression", "Interaction")),
  label = sprintf("p-value = %.3g\nrho = %.2f", trend_df$kw_p_value, trend_df$spearman_rho)
)

response_colors <- c(CR = "#4393c3", PR = "#8fb3d1", SD = "#dfa092", PD = "#d6604d")

p <- ggplot(plot_df, aes(x = Response, y = score, fill = Response)) +
  geom_boxplot(width = 0.65, alpha = 0.78, outlier.shape = NA, linewidth = 0.25) +
  geom_jitter(width = 0.14, size = 0.3, stroke = 0, alpha = 0.5) +
  geom_text(data = label_df, aes(x = 2.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.1, size = 5 / .pt, lineheight = 0.9) +
  facet_grid(score_type ~ module_label, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  scale_fill_manual(values = response_colors) +
  labs(
    title = "Response associations of predefined immune modules at the gene and interaction levels",
    x = NULL,
    y = "Module score"
  ) +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5, face = "plain", colour = "black"),
    plot.title = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, colour = "black"),
    strip.text = element_text(size = 5, colour = "black", margin = margin(1, 1, 1, 1)),
    strip.background = element_rect(fill = "grey95", linewidth = 0.3),
    panel.border = element_rect(linewidth = 0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(0.8, "mm"),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "figure7a.pdf"), p, width = 16, height = 8, units = "cm")

print(trend_df)
