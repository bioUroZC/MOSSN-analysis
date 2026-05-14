# ==============================================================================
# 1. Clear and load
# ==============================================================================
rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(data.table)

# ==============================================================================
# 2. Load data
# ==============================================================================

dieasename <- "Paclitaxel"
dataset_name <- "GSE66305"


base_path <- file.path("/proj/c.zihao/work1/2drugs/",
                       dieasename, "/", 
                       dataset_name, "data")

save_path <- file.path("/proj/c.zihao/work1/2drugs/",
                       dieasename, "/", 
                       dataset_name, "CSN")

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

setwd(base_path)
expr <- read.csv("exprSet_filtered.csv", header = TRUE, row.names = 1)
interGenes <- rownames(expr)

# Edge list
edge_raw <- read.csv("/proj/c.zihao/work1/1survival/links.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
edge_list <- edge_raw %>%
  dplyr::select(protein1, protein2) %>%
  dplyr::filter(protein1 %in% interGenes, protein2 %in% interGenes) %>%
  dplyr::distinct()
colnames(edge_list) <- c("Gene1", "Gene2")
rownames(edge_list) <- NULL
cat("Filtered edge list:", nrow(edge_list), "edges\n")

# ==============================================================================
# 3. Modular functions (already defined)
# ==============================================================================

compute_gene_boundaries <- function(data, boxsize = 0.1) {
  n_genes <- nrow(data)
  n_samples <- ncol(data)
  upper <- matrix(0, n_genes, n_samples)
  lower <- matrix(0, n_genes, n_samples)
  for (i in 1:n_genes) {
    g <- as.numeric(data[i, ])
    s <- sort(g, index.return = TRUE)
    s1 <- s$x
    s2 <- s$ix
    n0 <- n_samples - sum(sign(s1))
    h <- round(boxsize / 2 * sum(sign(s1)))
    k <- 1
    while (k <= n_samples) {
      s <- 0
      while ((k + s + 1) <= n_samples && s1[k + s + 1] == s1[k]) {
        s <- s + 1
      }
      if (s >= h) {
        upper[i, s2[k:(k + s)]] <- g[s2[k]]
        lower[i, s2[k:(k + s)]] <- g[s2[k]]
      } else {
        upper[i, s2[k:(k + s)]] <- g[s2[min(n_samples, k + s + h)]]
        lower[i, s2[k:(k + s)]] <- g[s2[max(ifelse(n0 > h, n0, 1), k - h)]]
      }
      k <- k + s + 1
    }
  }
  return(list(upper = upper, lower = lower))
}

compute_binary_matrix <- function(data, sample_index, upper, lower) {
  n_genes <- nrow(data)
  n_samples <- ncol(data)
  B <- matrix(0, n_genes, n_samples)
  for (j in 1:n_samples) {
    B[, j] <- (data[, j] <= upper[, sample_index]) & 
      (data[, j] >= lower[, sample_index])
  }
  return(B)
}

compute_edge_status <- function(B, edge_list, gene_index, alpha = 0.01) {
  a <- rowSums(B)
  n <- ncol(B)
  threshold <- -qnorm(alpha)
  edge_status <- integer(nrow(edge_list))
  for (i in 1:nrow(edge_list)) {
    g1 <- edge_list$Gene1[i]
    g2 <- edge_list$Gene2[i]
    idx1 <- gene_index[[g1]]
    idx2 <- gene_index[[g2]]
    if (is.null(idx1) || is.null(idx2)) {
      edge_status[i] <- NA
      next
    }
    score <- (sum(B[idx1, ] & B[idx2, ]) * n - a[idx1] * a[idx2]) /
      sqrt(a[idx1] * a[idx2] * (n - a[idx1]) * (n - a[idx2]) / (n - 1) + 1e-6)
    edge_status[i] <- as.integer(score > threshold)
  }
  return(edge_status)
}

# ==============================================================================
# 4. Compute CSN edges for all samples using your modular functions
# ==============================================================================
gene_index <- setNames(1:nrow(expr), rownames(expr))
sample_names <- colnames(expr)
edge_names <- paste(edge_list$Gene1, edge_list$Gene2, sep = "_")

# Step 1: boundaries (common for all samples)
cat("Computing gene boundaries...\n")
bounds <- compute_gene_boundaries(expr, boxsize = 0.1)

# Step 2: loop through each sample
result_list <- list()
for (i in seq_along(sample_names)) {
  cat("Processing sample", i, ":", sample_names[i], "\n")
  B <- compute_binary_matrix(expr, i, bounds$upper, bounds$lower)
  edge_status <- compute_edge_status(B, edge_list, gene_index, alpha = 0.01)
  result_list[[i]] <- edge_status
}

# Combine into matrix
result_matrix <- do.call(cbind, result_list)
colnames(result_matrix) <- sample_names
rownames(result_matrix) <- edge_names

# Convert to dataframe
result_df <- data.frame(edge = edge_names, result_matrix, row.names = NULL)

# ==============================================================================
# 5. Output
# ==============================================================================

setwd(save_path)

write.csv(result_df, "result_matrix.csv", row.names = FALSE) 

cat("All samples for CSN completed")

