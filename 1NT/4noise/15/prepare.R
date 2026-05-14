rm(list = ls())

set.seed(123)

base_dir <- "/proj/c.zihao/work1/1NT/4noise/15"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

noise_level <- 0.15

expr_file <- "/proj/c.zihao/work1/1NT/1data/exprset/LUAD_exprSet_filtered.csv"
expr <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)

sample_info <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
sample_info$group <- substr(sample_info$sample, 14, 16)
tumor_samples <- sample_info$sample[sample_info$group == "01A"]

expr_tumor <- expr[, tumor_samples, drop = FALSE]

gene_sd <- apply(expr_tumor, 1, sd, na.rm = TRUE)
noise_mat <- matrix(
  rnorm(length(expr_tumor), mean = 0, sd = rep(noise_level * gene_sd, ncol(expr_tumor))),
  nrow = nrow(expr_tumor),
  ncol = ncol(expr_tumor),
  byrow = FALSE
)
rownames(noise_mat) <- rownames(expr_tumor)
colnames(noise_mat) <- colnames(expr_tumor)

expr_noise <- expr_tumor + noise_mat
expr_noise[expr_noise < 0] <- 0

write.csv(
  expr_noise,
  file = file.path(data_dir, "LUAD_exprSet_noise.csv"),
  quote = FALSE
)

write.csv(
  data.frame(
    sample = tumor_samples,
    group = "01A",
    stringsAsFactors = FALSE
  ),
  file = file.path(data_dir, "LUAD_noise_samples.csv"),
  row.names = FALSE,
  quote = FALSE
)

normal_file <- "/proj/c.zihao/work1/0ref/Test8/combined_expr_df.csv"
normal_df <- read.csv(normal_file, header = TRUE, check.names = FALSE)
lung_normal_df <- normal_df[normal_df$organ == "Lung", , drop = FALSE]

write.csv(
  lung_normal_df,
  file = file.path(data_dir, "Lung_normal_reference.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    sample = lung_normal_df[[1]],
    organ = lung_normal_df$organ,
    stringsAsFactors = FALSE
  ),
  file = file.path(data_dir, "Lung_normal_reference_samples.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Noise level:", noise_level, "\n")
cat("Original samples:", ncol(expr), "\n")
cat("Tumor samples used:", ncol(expr_tumor), "\n")
cat("Saved noisy tumor expression matrix to:", file.path(data_dir, "LUAD_exprSet_noise.csv"), "\n")
cat("Saved unchanged Lung normal reference to:", file.path(data_dir, "Lung_normal_reference.csv"), "\n")
