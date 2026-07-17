rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

base_dir <- paste0(PROJ_ROOT, "/3drugs")
out_dir <- file.path(base_dir, "immune/IM210/case/results")

cat("Step 1: load score table\n")
score_wide <- read.csv(
  file.path(out_dir, "module_scores_wide.csv"),
  stringsAsFactors = FALSE
)

surv_df <- score_wide %>%
  filter(!is.na(Time), !is.na(OS)) %>%
  mutate(
    Time = as.numeric(Time),
    OS = as.numeric(OS)
  ) %>%
  filter(
    is.finite(Time),
    is.finite(OS),
    Time > 0
  )

cat("Samples with survival data:", nrow(surv_df), "\n")
cat("Events:", sum(surv_df$OS == 1), "\n")

cat("Step 2: prepare Kaplan-Meier groups\n")
surv_df <- surv_df %>%
  mutate(
    checkpoint_interaction_group = ifelse(
      int_Checkpoint_neighborhood >= median(int_Checkpoint_neighborhood, na.rm = TRUE),
      "High",
      "Low"
    ),
    checkpoint_expression_group = ifelse(
      expr_Checkpoint_neighborhood >= median(expr_Checkpoint_neighborhood, na.rm = TRUE),
      "High",
      "Low"
    )
  )

surv_df$checkpoint_interaction_group <- factor(surv_df$checkpoint_interaction_group, levels = c("Low", "High"))
surv_df$checkpoint_expression_group <- factor(surv_df$checkpoint_expression_group, levels = c("Low", "High"))

cat("Step 4: fit Kaplan-Meier curves\n")
fit_checkpoint_int <- survfit(
  Surv(Time, OS) ~ checkpoint_interaction_group,
  data = surv_df
)
fit_checkpoint_expr <- survfit(
  Surv(Time, OS) ~ checkpoint_expression_group,
  data = surv_df
)

p_checkpoint_int <- ggsurvplot(
  fit_checkpoint_int,
  data = surv_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#FDAE61", "#D7191C"),
  title = "Overall survival by checkpoint interaction score",
  legend.title = NULL,
  legend.labs = c("Low", "High"),
  xlab = "Time",
  ylab = "Overall survival probability"
)

p_checkpoint_expr <- ggsurvplot(
  fit_checkpoint_expr,
  data = surv_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#FDAE61", "#E6550D"),
  title = "Overall survival by checkpoint expression score",
  legend.title = NULL,
  legend.labs = c("Low", "High"),
  xlab = "Time",
  ylab = "Overall survival probability"
)

pdf(file.path(out_dir, "km_interaction.pdf"), width = 6.8, height = 7)
print(p_checkpoint_int)
dev.off()

pdf(file.path(out_dir, "km_expression.pdf"), width = 6.8, height = 7)
print(p_checkpoint_expr)
dev.off()
cat("Finished survival step.\n")
