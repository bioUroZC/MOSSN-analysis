rm(list = ls())

library(dplyr)
library(ggplot2)

base_dir <- "/proj/c.zihao/work1/3drugs"
out_dir <- file.path(base_dir, "immune/IM210/case/results")

module_order <- c(
  "IFNG_response",
  "Antigen_presentation",
  "Cytotoxicity",
  "Checkpoint_neighborhood"
)

pretty_module <- c(
  IFNG_response = "IFN-gamma response",
  Antigen_presentation = "Antigen presentation",
  Cytotoxicity = "Cytotoxicity",
  Checkpoint_neighborhood = "Checkpoint neighborhood"
)

cat("Step 1: load module scores\n")
score_long <- read.csv(
  file.path(out_dir, "module_scores_long.csv"),
  stringsAsFactors = FALSE
)
score_wide <- read.csv(
  file.path(out_dir, "module_scores_wide.csv"),
  stringsAsFactors = FALSE
)

meta_df <- score_wide[, c("Sample", "Response", "response_ordinal")]

score_long <- score_long %>%
  left_join(meta_df, by = "Sample") %>%
  filter(
    module %in% module_order,
    Response %in% c("CR", "PR", "SD", "PD")
  )

score_long$Response <- factor(score_long$Response, levels = c("CR", "PR", "SD", "PD"))

cat("Step 2: run Kruskal-Wallis and Spearman trend tests\n")
trend_df <- bind_rows(
  lapply(module_order, function(module_name) {
    sub_df <- score_long %>%
      filter(module == module_name)

    kw_expr <- kruskal.test(expression_score ~ Response, data = sub_df)
    kw_int <- kruskal.test(interaction_score ~ Response, data = sub_df)

    sp_expr <- suppressWarnings(cor.test(
      sub_df$expression_score,
      sub_df$response_ordinal,
      method = "spearman",
      exact = FALSE
    ))
    sp_int <- suppressWarnings(cor.test(
      sub_df$interaction_score,
      sub_df$response_ordinal,
      method = "spearman",
      exact = FALSE
    ))

    data.frame(
      module = module_name,
      score_type = c("Expression", "Interaction"),
      n = nrow(sub_df),
      median_CR = c(
        median(sub_df$expression_score[sub_df$Response == "CR"], na.rm = TRUE),
        median(sub_df$interaction_score[sub_df$Response == "CR"], na.rm = TRUE)
      ),
      median_PR = c(
        median(sub_df$expression_score[sub_df$Response == "PR"], na.rm = TRUE),
        median(sub_df$interaction_score[sub_df$Response == "PR"], na.rm = TRUE)
      ),
      median_SD = c(
        median(sub_df$expression_score[sub_df$Response == "SD"], na.rm = TRUE),
        median(sub_df$interaction_score[sub_df$Response == "SD"], na.rm = TRUE)
      ),
      median_PD = c(
        median(sub_df$expression_score[sub_df$Response == "PD"], na.rm = TRUE),
        median(sub_df$interaction_score[sub_df$Response == "PD"], na.rm = TRUE)
      ),
      kw_p_value = c(kw_expr$p.value, kw_int$p.value),
      spearman_rho = c(as.numeric(sp_expr$estimate), as.numeric(sp_int$estimate)),
      spearman_p_value = c(sp_expr$p.value, sp_int$p.value),
      stringsAsFactors = FALSE
    )
  })
)

trend_df$kw_fdr <- p.adjust(trend_df$kw_p_value, method = "BH")
trend_df$spearman_fdr <- p.adjust(trend_df$spearman_p_value, method = "BH")
trend_df$module_label <- pretty_module[trend_df$module]
trend_df <- trend_df %>%
  arrange(match(module, module_order), match(score_type, c("Expression", "Interaction")))

cat("Step 3: prepare plotting tables\n")
plot_df <- bind_rows(
  score_long %>%
    transmute(
      Sample = Sample,
      module = module,
      module_label = pretty_module[module],
      Response = Response,
      score_type = "Expression",
      score = expression_score
    ),
  score_long %>%
    transmute(
      Sample = Sample,
      module = module,
      module_label = pretty_module[module],
      Response = Response,
      score_type = "Interaction",
      score = interaction_score
    )
)

plot_df$module_label <- factor(plot_df$module_label, levels = pretty_module[module_order])
plot_df$score_type <- factor(plot_df$score_type, levels = c("Expression", "Interaction"))
plot_df$Response <- factor(plot_df$Response, levels = c("CR", "PR", "SD", "PD"))

label_df <- trend_df %>%
  transmute(
    module_label = pretty_module[module],
    score_type = score_type,
    label = sprintf("KW p = %.3g\nrho = %.2f", kw_p_value, spearman_rho)
  )
label_df$module_label <- factor(label_df$module_label, levels = pretty_module[module_order])
label_df$score_type <- factor(label_df$score_type, levels = c("Expression", "Interaction"))

cat("Step 4: draw trend plots\n")
p_trend <- ggplot(
  plot_df,
  aes(x = Response, y = score, fill = Response)
) +
  geom_boxplot(width = 0.65, alpha = 0.78, outlier.shape = NA) +
  geom_jitter(width = 0.14, size = 0.9, alpha = 0.45) +
  facet_grid(score_type ~ module_label, scales = "free_y") +
  scale_fill_manual(values = c(CR = "#2C7BB6", PR = "#ABD9E9", SD = "#FDAE61", PD = "#D7191C")) +
  labs(
    title = "Predefined immune module scores track the clinical response gradient",
    subtitle = "CR to PD ordering was assessed without response-driven feature selection",
    x = NULL,
    y = "Module score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  ) +
  geom_text(
    data = label_df,
    aes(x = 2.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.1,
    size = 3
  )

ggsave(
  file.path(out_dir, "response_trend_boxplots.pdf"),
  p_trend,
  width = 11,
  height = 6
)
cat("Finished response-trend step.\n")
print(trend_df[, c("module", "score_type", "kw_p_value", "kw_fdr", "spearman_rho", "spearman_p_value", "spearman_fdr")])
