rm(list = ls())

library(ggplot2)
library(data.table)

method_dir    <- "/proj/c.zihao/work1/1NT/2string"
metadata_path <- "/proj/c.zihao/work1/1NT/1data/TCGA/metadata.csv"
out_dir       <- "/proj/c.zihao/work1/1NT/2string/1distance/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method    <- "MOSSN_uniform"
cancer    <- "LUAD"
top_frac  <- 0.20
seed      <- 1

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

mat <- fread(file.path(method_dir, method, "merged_matrix.csv")) |> as.data.frame()
id_col <- intersect(c("Interaction", "V1", "Unnamed: 0"), colnames(mat))[1]
rownames(mat) <- mat[[id_col]]
mat[[id_col]] <- NULL
mat <- abs(mat)

cols    <- intersect(metadata$Sample[metadata$Type == cancer], colnames(mat))
sub_mat <- mat[, cols, drop = FALSE]
sub_mat[is.na(sub_mat)] <- 0

cv    <- apply(sub_mat, 1, sd) / abs(rowMeans(sub_mat))
cv    <- cv[is.finite(cv)]
keep  <- names(sort(cv, decreasing = TRUE))[seq_len(floor(length(cv) * top_frac))]

X_raw <- t(as.matrix(sub_mat[keep, , drop = FALSE]))
X_raw <- X_raw[, apply(X_raw, 2, sd) > 0, drop = FALSE]
X     <- scale(X_raw)

set.seed(seed)
pca     <- prcomp(X, center = FALSE, scale. = FALSE)
var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100

df <- data.frame(
  PC1  = pca$x[, 1],
  PC2  = pca$x[, 2],
  type = ifelse(grepl("11A$", rownames(X)), "Normal", "Tumor")
)

p <- ggplot(df, aes(PC1, PC2, color = type)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(type = "norm", linewidth = 0.6, linetype = "dashed", level = 0.90) +
  scale_color_manual(values = c(Tumor = "#d62728", Normal = "#1f77b4")) +
  labs(
    title = "LUAD - MOSSN",
    x     = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y     = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position  = "bottom"
  )

ggsave(file.path(out_dir, "pca_LUAD_MOSSN.pdf"), p,
       width = 6, height = 6)
