LIONcal <- function (exprSetFile, ppiFile){ 

library(dplyr)
exp <- read.csv(exprSetFile, header = TRUE, row.names = 1)
exp <- exp[apply(exp, 1, sd) > 0, ]

link <- read.csv(ppiFile, header = TRUE, row.names = 1)
print(dim(link))

# ===================================================

# ===================================================

common_genes <- intersect(rownames(exp), unique(c(link$protein1, link$protein2)))

link_filtered <- link %>%
  dplyr::filter(protein1 %in% common_genes & protein2 %in% common_genes) %>%
  dplyr::distinct(protein1, protein2, .keep_all = TRUE)


used_genes <- unique(c(link_filtered$protein1, link_filtered$protein2))
dat <- exp[rownames(exp) %in% used_genes, ]
dat <- as.matrix(dat)

cat("Number of genes used: ", nrow(dat), "\n")
cat("Number of edges retained: ", nrow(link_filtered), "\n")


# Custom network function returning only selected edge correlations
netFun_string <- function(x) {
  x <- x[apply(x, 1, sd) > 0, ]
  if (nrow(x) < 2) return(rep(NA, nrow(link_filtered)))
  
  cor_mat <- cor(t(x), method = "pearson")
  
  edge_vals <- mapply(function(g1, g2) {
    if (g1 %in% rownames(x) && g2 %in% rownames(x)) {
      return(cor_mat[g1, g2])
    } else {
      return(NA)
    }
  }, link_filtered$protein1, link_filtered$protein2)
  
  names(edge_vals) <- paste0(link_filtered$protein1, "_", link_filtered$protein2)
  return(edge_vals)
}

# Lightweight LIONESS implementation with progress messages
lioness_vector <- function(x, f = netFun_string, ...) {
  if (!is.function(f)) stop("Function required")
  if (!is.matrix(x)) stop("Input must be a numeric matrix")
  
  nrsamples <- ncol(x)
  samples <- colnames(x)
  
  cat("Computing aggregate network...\n")
  net_agg <- f(x)
  edge_names <- names(net_agg)
  if (is.null(edge_names)) stop("netFun output must be a named vector")
  
  lioness_mat <- matrix(NA, nrow = length(net_agg), ncol = nrsamples)
  colnames(lioness_mat) <- samples
  rownames(lioness_mat) <- edge_names
  
  cat("Running LIONESS leave-one-out...\n")
  for (i in 1:nrsamples) {
    cat(sprintf("Processing sample %d (%s) of %d\n", i, samples[i], nrsamples))
    x_minus_i <- x[, -i, drop = FALSE]
    net_minus_i <- f(x_minus_i)
    lioness_mat[, i] <- nrsamples * (net_agg - net_minus_i) + net_minus_i
  }
  
  cat("LIONESS complete. Formatting output...\n")
  
  edge_df <- data.frame(edge = rownames(lioness_mat), lioness_mat, check.names = FALSE)
  edge_df <- edge_df %>%
    tidyr::separate(edge, into = c("reg", "tar"), sep = "_") %>%
    dplyr::select(reg, tar, everything())
  
  return(edge_df)
}

# Run LIONESS
result <- lioness_vector(dat, netFun_string)

return(result)

}
