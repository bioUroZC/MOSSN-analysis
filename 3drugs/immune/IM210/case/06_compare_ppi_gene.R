rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(ggplot2)
library(pROC)

base_dir <- paste0(PROJ_ROOT, "/3drugs")
out_dir <- file.path(base_dir, "immune/IM210/case/results")
cat("Step 1: load follow-up data\n")

follow_up <- read.csv(
  file.path(base_dir, "immune/IM210/data/IMvigor210_FollowUp.csv"),
  header = TRUE,
  row.names = 1
)
follow_up$Sample <- rownames(follow_up)
follow_up <- subset(
  follow_up,
  select = c(
    "Sample",
    "Best.Confirmed.Overall.Response",
    "binaryResponse",
    "Immune.phenotype",
    "Tissue",
    "os",
    "censOS"
  )
)
names(follow_up) <- c("Sample", "Response", "Binary", "Phenotype", "Tissue", "Time", "OS")
cat("Step 2: load gene expression\n")

gene_expr <- read.csv(
  file.path(base_dir, "immune/IM210/data/exprSet_filtered.csv"),
  header = TRUE,
  row.names = 1
)

gene_expr <- apply(gene_expr, 2, function(x) {
  if (length(unique(x)) == 1) {
    return(rep(0, length(x)))
  }
  ranks <- rank(x, ties.method = "average")
  (ranks - 1) / (length(ranks) - 1)
})

gene_expr <- gene_expr[rownames(gene_expr) %in% c("PDCD1", "CD274"), , drop = FALSE]
gene_expr <- as.data.frame(t(gene_expr))
gene_expr$Sample <- rownames(gene_expr)
gene_expr$Gene_max <- pmax(gene_expr$PDCD1, gene_expr$CD274, na.rm = TRUE)
gene_expr$Gene_mean <- rowMeans(gene_expr[, c("PDCD1", "CD274")], na.rm = TRUE)
cat("Step 3: load link matrix\n")

link_path <- file.path(base_dir, "immune/1Matrix/MOSSN_uniform_Matrix.csv")
con <- file(link_path, open = "r")
header_line <- readLines(con, n = 1)
link_lines <- character(0)
repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  if (grepl('^"?(CD274_PDCD1|CD274_CD8A)"?,', line)) {
    link_lines <- c(link_lines, line)
  }
  if (length(link_lines) == 2) {
    break
  }
}
close(con)
link_expr <- read.csv(
  text = c(header_line, link_lines),
  header = TRUE,
  check.names = FALSE
)
rownames(link_expr) <- link_expr$Interaction
link_expr$Interaction <- NULL
link_expr <- as.data.frame(t(link_expr))
link_expr$Sample <- rownames(link_expr)
cat("Step 4: build matched cohort\n")

dat <- follow_up %>%
  inner_join(gene_expr, by = "Sample") %>%
  inner_join(link_expr, by = "Sample") %>%
  filter(
    Tissue == "bladder",
    Response %in% c("CR", "PD")
  ) %>%
  mutate(responder = ifelse(Response == "CR", 1, 0)) %>%
  select(
    Sample,
    Response,
    responder,
    PDCD1,
    CD274,
    Gene_max,
    Gene_mean,
    CD274_PDCD1,
    CD274_CD8A
  ) %>%
  filter(complete.cases(.))
cat("Matched samples:", nrow(dat), "\n")

score_map <- c(
  PDCD1 = "PDCD1",
  CD274 = "CD274",
  Gene_max = "Gene_max",
  Gene_mean = "Gene_mean",
  PPI_link = "CD274_PDCD1"
)

roc_list <- lapply(score_map, function(col_name) {
  roc(response = dat$responder, predictor = dat[[col_name]], quiet = TRUE)
})
cat("Step 5: ROC objects done\n")

summary_df <- do.call(
  rbind,
  lapply(names(score_map), function(model_name) {
    roc_obj <- roc_list[[model_name]]
    ci_obj <- ci.auc(roc_obj, method = "delong")
    data.frame(
      model = model_name,
      feature = score_map[[model_name]],
      n = nrow(dat),
      responders = sum(dat$responder == 1),
      nonresponders = sum(dat$responder == 0),
      auc = as.numeric(auc(roc_obj)),
      ci_low = as.numeric(ci_obj[1]),
      ci_high = as.numeric(ci_obj[3]),
      stringsAsFactors = FALSE
    )
  })
)
cat("Step 6: summary table done\n")

delong_df <- data.frame(
  comparison = c("PPI_link_vs_Gene_max", "PPI_link_vs_PDCD1", "PPI_link_vs_CD274"),
  p_value = c(
    roc.test(roc_list$Gene_max, roc_list$PPI_link, method = "delong", alternative = "less")$p.value,
    roc.test(roc_list$PDCD1, roc_list$PPI_link, method = "delong", alternative = "less")$p.value,
    roc.test(roc_list$CD274, roc_list$PPI_link, method = "delong", alternative = "less")$p.value
  ),
  stringsAsFactors = FALSE
)
cat("Step 7: DeLong tests done\n")

cat("Step 8: build plot\n")

make_label <- function(model_name, short = FALSE) {
  row <- summary_df[summary_df$model == model_name, ]
  if (short) {
    sprintf("%s (AUC %.3f)", model_name, row$auc)
  } else {
    sprintf(
      "%s (AUC %.3f, 95%% CI %.3f-%.3f)",
      model_name, row$auc, row$ci_low, row$ci_high
    )
  }
}

p <- ggroc(
  list(
    `Gene max` = roc_list$Gene_max,
    `PPI link` = roc_list$PPI_link
  ),
  legacy.axes = TRUE,
  size = 1.2
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey60"
  ) +
  coord_equal() +
  labs(
    title = "CD274 & PDCD1: interaction edge vs gene expression",
    subtitle = sprintf(
      "DeLong p = %.4g (CD274-PDCD1 edge vs max expression)",
      delong_df$p_value[delong_df$comparison == "PPI_link_vs_Gene_max"]
    ),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  scale_color_manual(
    values = c(
      `Gene max` = "#2171B5",
      `PPI link` = "#CB181D"
    ),
    labels = c(
      make_label("Gene_max", short = TRUE),
      make_label("PPI_link", short = TRUE)
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.64, 0.20),
    legend.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "ppi_vs_gene_roc.pdf"),
  p,
  width = 7,
  height = 5.5
)
cat("Step 9: figure written\n")

print(summary_df)
print(delong_df)
cat("Finished. Results written to:", out_dir, "\n")
