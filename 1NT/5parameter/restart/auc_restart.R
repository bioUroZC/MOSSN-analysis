rm(list = ls())

library(data.table)
library(dplyr)

base_dir <- "/proj/c.zihao/work1/1NT/5parameter"
restart_dir <- file.path(base_dir, "restart")
out_dir <- file.path(restart_dir, "auc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_frac <- 0.20
nstart <- 50
seed <- 123

metadata <- read.csv(file.path(base_dir, "data", "LUAD_paired_metadata.csv"), stringsAsFactors = FALSE)
sample_ids <- metadata$sample

calc_auc <- function(labels, scores) {
  labels <- as.integer(labels)
  pos <- sum(labels == 1)
  neg <- sum(labels == 0)
  if (pos == 0 || neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1]) - pos * (pos + 1) / 2) / (pos * neg)
}

calc_aupr <- function(labels, scores) {
  labels <- as.integer(labels)
  pos <- sum(labels == 1)
  if (pos == 0) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  y <- labels[ord]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  precision <- tp / (tp + fp)
  recall <- tp / pos
  precision_prev <- c(1, precision[-length(precision)])
  recall_prev <- c(0, recall[-length(recall)])
  sum((recall - recall_prev) * (precision + precision_prev) / 2)
}

safe_div <- function(num, den) {
  ifelse(den == 0, NA_real_, num / den)
}

read_matrix_file <- function(matrix_file) {
  pair_matrix <- fread(matrix_file) |> as.data.frame()
  id_col <- intersect(c("Interaction", "V1", "Unnamed: 0"), colnames(pair_matrix))
  if (length(id_col) == 0) {
    stop("Unable to find interaction ID column in: ", matrix_file)
  }
  rownames(pair_matrix) <- pair_matrix[[id_col[1]]]
  pair_matrix[[id_col[1]]] <- NULL
  pair_matrix
}

analyze_matrix <- function(matrix_file, sample_ids, top_frac, nstart, seed) {
  pair_matrix <- read_matrix_file(matrix_file)

  common_samples <- intersect(sample_ids, colnames(pair_matrix))
  if (length(common_samples) < 4) {
    return(NULL)
  }

  data_final <- pair_matrix[, common_samples, drop = FALSE]
  data_final[is.na(data_final)] <- 0
  data_final <- abs(data_final)

  means <- rowMeans(data_final)
  sds <- apply(data_final, 1, sd)
  cv <- sds / ifelse(means == 0, NA, abs(means))
  cv <- stats::na.omit(cv)
  if (length(cv) == 0) {
    return(NULL)
  }

  top_n <- max(1, ceiling(top_frac * length(cv)))
  keep <- names(sort(cv, decreasing = TRUE))[seq_len(min(top_n, length(cv)))]
  data_top_n <- data_final[keep, , drop = FALSE]

  raw_matrix <- t(as.matrix(data_top_n))
  feature_sd <- apply(raw_matrix, 2, sd)
  raw_matrix <- raw_matrix[, feature_sd > 0, drop = FALSE]
  if (ncol(raw_matrix) == 0) {
    return(NULL)
  }

  X <- scale(raw_matrix)
  sample_groups <- ifelse(grepl("11A$", rownames(X)), "Normal", "Tumor")
  label_pos <- as.integer(sample_groups == "Tumor")

  set.seed(seed)
  km <- kmeans(X, centers = 2, nstart = nstart)
  pred_if_c1 <- ifelse(km$cluster == 1, "Normal", "Tumor")
  pred_if_c2 <- ifelse(km$cluster == 2, "Normal", "Tumor")
  acc1 <- mean(pred_if_c1 == sample_groups)
  acc2 <- mean(pred_if_c2 == sample_groups)
  pred_groups <- if (acc1 >= acc2) pred_if_c1 else pred_if_c2
  pred_pos <- as.integer(pred_groups == "Tumor")

  accuracy <- mean(pred_groups == sample_groups)
  tp <- sum(pred_pos == 1 & label_pos == 1)
  fp <- sum(pred_pos == 1 & label_pos == 0)
  fn <- sum(pred_pos == 0 & label_pos == 1)
  precision <- safe_div(tp, tp + fp)
  recall <- safe_div(tp, tp + fn)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }

  pca_res <- prcomp(X, center = FALSE, scale. = FALSE)
  pc1 <- pca_res$x[, 1]
  if (mean(pc1[label_pos == 1]) < mean(pc1[label_pos == 0])) {
    pc1 <- -pc1
  }
  auc <- calc_auc(label_pos, pc1)
  aupr <- calc_aupr(label_pos, pc1)

  data.frame(
    n_sample = nrow(X),
    n_variable_features = length(cv),
    n_selected_features = length(keep),
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1 = f1,
    auc = auc,
    aupr = aupr,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()

alpha_dirs <- list.dirs(restart_dir, recursive = FALSE, full.names = TRUE)
alpha_dirs <- alpha_dirs[basename(alpha_dirs) != "auc"]
alpha_dirs <- alpha_dirs[basename(alpha_dirs) != "distance"]
alpha_dirs <- alpha_dirs[basename(alpha_dirs) != "plots"]

for (alpha_dir in alpha_dirs) {
  alpha_name <- basename(alpha_dir)
  matrix_file <- file.path(alpha_dir, "merged_matrix.csv")
  if (!file.exists(matrix_file)) next

  summary_row <- analyze_matrix(matrix_file, sample_ids, top_frac, nstart, seed)
  if (is.null(summary_row)) next

  summary_row$alpha <- alpha_name
  summary_rows[[alpha_name]] <- summary_row

  cat(
    alpha_name,
    ": n_sample=", summary_row$n_sample,
    " selected=", summary_row$n_selected_features,
    " accuracy=", round(summary_row$accuracy, 3),
    " precision=", round(summary_row$precision, 3),
    " recall=", round(summary_row$recall, 3),
    " f1=", round(summary_row$f1, 3),
    " auc=", round(summary_row$auc, 3),
    " aupr=", round(summary_row$aupr, 3),
    "\n"
  )
}

summary_df <- bind_rows(summary_rows) |>
  select(
    alpha, n_sample, n_variable_features, n_selected_features,
    accuracy, precision, recall, f1, auc, aupr
  )

write.csv(summary_df, file.path(out_dir, "restart_metrics_summary.csv"), row.names = FALSE)
