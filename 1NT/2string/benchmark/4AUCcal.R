rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(data.table)

base_dir <- paste0(PROJ_ROOT, "/1NT/2string/benchmark")
metadata_path <- paste0(PROJ_ROOT, "/1NT/1data/TCGA/metadata.csv")
output_dir <- file.path(base_dir, "1result")
dir.create(output_dir, showWarnings = FALSE)

top_fracs <- c(0.05, 0.10, 0.15, 0.20)
k <- 2
nstart <- 50
seed <- 1

cancers <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
    "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

method_files <- sort(list.files(
  base_dir,
  pattern = "merged_matrix\\.csv$",
  recursive = TRUE,
  full.names = TRUE
))
method_names <- basename(dirname(method_files))

cat("Methods:\n")
cat(paste(" -", method_names, collapse = "\n"), "\n")

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

make_metric_summary <- function(df, metric_name, frac_label) {
  sub_df <- df[df$feature_fraction == frac_label, c("cancer", "method", metric_name)]
  out <- reshape(
    sub_df,
    idvar = "cancer",
    timevar = "method",
    direction = "wide"
  )
  colnames(out) <- gsub(paste0("^", metric_name, "\\."), "", colnames(out))
  out
}

cluster_df <- data.frame(
  method = character(),
  feature_fraction = character(),
  cancer = character(),
  n_sample = integer(),
  n_variable_features = integer(),
  n_selected_features = integer(),
  accuracy = numeric(),
  precision = numeric(),
  recall = numeric(),
  f1 = numeric(),
  auc = numeric(),
  aupr = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(method_files)) {
  cat("\n[Method]", method_names[i], "\n")

  mat <- fread(method_files[i]) |> as.data.frame()
  id_col <- intersect(c("Interaction", "V1", "Unnamed: 0"), colnames(mat))
  if (length(id_col) == 0) {
    stop("No interaction ID column found in: ", method_files[i])
  }
  rownames(mat) <- mat[[id_col[1]]]
  mat[[id_col[1]]] <- NULL
  mat <- abs(mat)

  for (cancer in cancers) {
    cancer_samples <- metadata$Sample[metadata$Type == cancer]
    cols <- intersect(cancer_samples, colnames(mat))
    if (length(cols) < 4) next

    sub_mat <- mat[, cols, drop = FALSE]
    sub_mat[is.na(sub_mat)] <- 0

    means <- rowMeans(sub_mat)
    sds <- apply(sub_mat, 1, sd)
    cv <- sds / ifelse(means == 0, NA, abs(means))
    cv <- na.omit(cv)
    if (length(cv) == 0) next

    for (top_frac in top_fracs) {
      frac_label <- paste0(as.integer(top_frac * 100), "%")
      top_n <- max(1, floor(length(cv) * top_frac))
      keep <- names(sort(cv, decreasing = TRUE))[seq_len(min(top_n, length(cv)))]
      sub_top <- sub_mat[keep, , drop = FALSE]

      X_raw <- t(as.matrix(sub_top))
      nonconst <- apply(X_raw, 2, sd) > 0
      X_raw <- X_raw[, nonconst, drop = FALSE]
      if (ncol(X_raw) == 0) {
        cat(" ", cancer, frac_label, ": skipped (all interactions constant)\n")
        next
      }
      X <- scale(X_raw)

      sample_ids <- rownames(X)
      sample_groups <- ifelse(grepl("11A$", sample_ids), "Normal", "Tumor")
      label_pos <- as.integer(sample_groups == "Tumor")

      set.seed(seed)
      km <- kmeans(X, centers = k, nstart = nstart)
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

      cluster_df <- rbind(cluster_df, data.frame(
        method = method_names[i],
        feature_fraction = frac_label,
        cancer = cancer,
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
      ))

      cat(
        " ", cancer, frac_label,
        ": selected=", length(keep),
        " accuracy=", round(accuracy, 3),
        " precision=", round(precision, 3),
        " recall=", round(recall, 3),
        " f1=", round(f1, 3),
        " auc=", round(auc, 3),
        " aupr=", round(aupr, 3),
        "\n"
      )
    }
  }
}

round_cols <- c("accuracy", "precision", "recall", "f1", "auc", "aupr")
cluster_df[round_cols] <- lapply(cluster_df[round_cols], round, digits = 5)

write.csv(
  cluster_df,
  file.path(output_dir, "results_Cluster.csv"),
  row.names = FALSE
)
cat("\nSaved to", file.path(output_dir, "results_Cluster.csv"), "\n")

metric_summaries <- lapply(setNames(paste0(as.integer(top_fracs * 100), "%"), paste0("top_", as.integer(top_fracs * 100))), function(frac_label) {
  list(
    accuracy = make_metric_summary(cluster_df, "accuracy", frac_label),
    precision = make_metric_summary(cluster_df, "precision", frac_label),
    recall = make_metric_summary(cluster_df, "recall", frac_label),
    f1 = make_metric_summary(cluster_df, "f1", frac_label),
    auc = make_metric_summary(cluster_df, "auc", frac_label),
    aupr = make_metric_summary(cluster_df, "aupr", frac_label)
  )
})
print(metric_summaries)

metric_cols <- c("accuracy", "precision", "recall", "f1", "auc", "aupr")
overall_metrics <- lapply(metric_cols, function(metric_name) {
  out <- aggregate(
    cluster_df[[metric_name]],
    by = list(method = cluster_df$method, feature_fraction = cluster_df$feature_fraction),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  colnames(out)[3] <- paste0("mean_", metric_name)
  out[order(out[[3]], decreasing = TRUE), ]
})
names(overall_metrics) <- metric_cols

print(overall_metrics)
