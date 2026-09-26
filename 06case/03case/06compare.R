rm(list = ls())

library(ggplot2)
library(pROC)

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")

clinical_file <- file.path(base_dir, "01data/IMvigor210_FollowUp.csv")
expr_file <- file.path(base_dir, "01data/exprSet_filtered.csv")
matrix_file <- file.path(base_dir, "04results/MOSSN/merged_matrix.csv")

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_rank <- apply(expr, 2, rank)
expr_rank <- (expr_rank - 1) / (nrow(expr) - 1)

gene_df <- data.frame(
  Sample = colnames(expr),
  PDCD1 = expr_rank["PDCD1", ],
  CD274 = expr_rank["CD274", ]
)
gene_df$Gene_max <- pmax(gene_df$PDCD1, gene_df$CD274)
gene_df$Gene_mean <- (gene_df$PDCD1 + gene_df$CD274) / 2

link <- data.table::fread(matrix_file, data.table = FALSE)
link <- link[link$Interaction == "CD274_PDCD1", ]
link_df <- data.frame(
  Sample = colnames(link)[-1],
  CD274_PDCD1 = as.numeric(link[1, -1])
)

clinical <- read.csv(clinical_file, row.names = 1)
dat <- data.frame(
  Sample = rownames(clinical),
  Response = clinical$Best.Confirmed.Overall.Response,
  Tissue = clinical$Tissue
)
dat <- dat[dat$Tissue %in% "bladder" & dat$Response %in% c("CR", "PD"), ]
dat <- dat[dat$Sample %in% gene_df$Sample & dat$Sample %in% link_df$Sample, ]
dat$responder <- ifelse(dat$Response == "CR", 1, 0)

dat <- cbind(dat, gene_df[match(dat$Sample, gene_df$Sample), -1])
dat$CD274_PDCD1 <- link_df$CD274_PDCD1[match(dat$Sample, link_df$Sample)]
cat("Matched samples:", nrow(dat), "\n")

features <- c(
  PDCD1 = "PDCD1",
  CD274 = "CD274",
  Gene_max = "Gene_max",
  Gene_mean = "Gene_mean",
  PPI_link = "CD274_PDCD1"
)

roc_list <- list()
summary_df <- data.frame()
for (model in names(features)) {
  roc_obj <- roc(response = dat$responder, predictor = dat[[features[model]]], quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "delong")
  roc_list[[model]] <- roc_obj

  summary_df <- rbind(summary_df, data.frame(
    model = model,
    feature = features[model],
    n = nrow(dat),
    responders = sum(dat$responder == 1),
    nonresponders = sum(dat$responder == 0),
    auc = as.numeric(auc(roc_obj)),
    ci_low = ci_obj[1],
    ci_high = ci_obj[3]
  ))
}
rownames(summary_df) <- NULL

delong_df <- data.frame()
for (model in c("Gene_max", "PDCD1", "CD274")) {
  test <- roc.test(roc_list[[model]], roc_list$PPI_link, method = "delong", alternative = "less")
  delong_df <- rbind(delong_df, data.frame(
    comparison = paste0("PPI_link_vs_", model),
    p_value = test$p.value
  ))
}

auc_gene <- summary_df$auc[summary_df$model == "Gene_max"]
auc_link <- summary_df$auc[summary_df$model == "PPI_link"]
p_value <- delong_df$p_value[delong_df$comparison == "PPI_link_vs_Gene_max"]

p <- ggroc(
  list(`Gene max` = roc_list$Gene_max, `PPI link` = roc_list$PPI_link),
  legacy.axes = TRUE,
  size = 0.4
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              linewidth = 0.3, color = "grey60") +
  coord_equal() +
  scale_color_manual(
    values = c(`Gene max` = "#4393c3", `PPI link` = "#d6604d"),
    labels = c(sprintf("Gene expression value (AUC %.3f)", auc_gene),
               sprintf("Interaction value (AUC %.3f)", auc_link))
  ) +
  labs(
    title = "ROC: gene vs interaction",
    subtitle = sprintf("p-value = %.4g", p_value),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw(base_size = 3) +
  theme(
    text = element_text(size = 3, face = "plain", colour = "black"),
    plot.title = element_text(size = 3),
    plot.subtitle = element_text(size = 3),
    axis.title = element_text(size = 3),
    axis.text = element_text(size = 3, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 3),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.02),
    legend.justification.inside = c(1, 0),
    legend.background = element_blank(),
    panel.border = element_rect(linewidth = 0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(0.8, "mm"),
    plot.margin = margin(1, 2, 1, 1)
  )

ggsave(file.path(out_dir, "figure7d.pdf"), p, width = 4, height = 4, units = "cm")

print(summary_df)
print(delong_df)
