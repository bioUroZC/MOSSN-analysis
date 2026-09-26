rm(list = ls())

library(ggplot2)
library(survival)
library(survminer)

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")

surv_df <- read.csv(file.path(out_dir, "module_scores_wide.csv"))
surv_df <- surv_df[!is.na(surv_df$Time) & !is.na(surv_df$OS) & surv_df$Time > 0, ]

cat("Samples with survival data:", nrow(surv_df), "\n")
cat("Events:", sum(surv_df$OS == 1), "\n")

figures <- data.frame(
  score_col = c("expr_Checkpoint_neighborhood", "int_Checkpoint_neighborhood"),
  title = c("Gene-based survival stratification", "Interaction-based survival stratification"),
  low_color = c("#d6604d", "#d6604d"),
  high_color = c("#4393c3", "#4393c3"),
  file = c("figure7b.pdf", "figure7c.pdf")
)

for (i in seq_len(nrow(figures))) {
  
  score <- surv_df[[figures$score_col[i]]]
  surv_df$group <- factor(ifelse(score >= median(score), "High", "Low"),
                          levels = c("Low", "High"))

  fit <- survfit(Surv(Time, OS) ~ group, data = surv_df)
  pval <- surv_pvalue(fit, data = surv_df)$pval  

  p <- ggsurvplot(
    fit,
    data = surv_df,
    risk.table = TRUE,
    risk.table.height = 0.3,
    risk.table.fontsize = 3 / .pt,
    pval = sprintf("p-value = %.3g", pval),
    pval.size = 3 / .pt,
    conf.int = FALSE,
    palette = c(figures$low_color[i], figures$high_color[i]),
    title = figures$title[i],
    legend.title = NULL,
    legend.labs = c("Low", "High"),
    xlab = "Time (months)",
    xlim = c(0, 20),
    break.time.by = 5,
    ylab = "Overall survival probability",
    size = 0.3,
    censor.size = 1.2
  )

  p$plot <- p$plot +
    theme(
      plot.title = element_text(size = 3, colour = "black"),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 3, colour = "black"),
      axis.title.y = element_text(size = 3, colour = "black"),
      axis.text.x = element_text(size = 3, colour = "black"),
      axis.text.y = element_text(size = 3, colour = "black"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(0.8, "mm"),
      legend.text = element_text(size = 3, colour = "black"),
      legend.title = element_blank(),
      legend.key.size = unit(2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-4, 0, -6, 0),
      plot.margin = margin(1, 2, 0, 1)
    )

  p$table <- p$table +
    theme(
      plot.title = element_text(size = 3, colour = "black"),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 3, colour = "black"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 3, colour = "black"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(0.8, "mm"),
      plot.margin = margin(0, 2, 1, 1)
    )
  
  p$table$theme$axis.text.y$size <- 3

  pdf(file.path(out_dir, figures$file[i]), width = 4 / 2.54, height = 4 / 2.54)
  print(p, newpage = FALSE)
  dev.off()
}
