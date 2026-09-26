rm(list = ls())

base_dir <- "/proj/c.zihao/work1/02analysis/02noise"
LEVELS <- c("05", "10", "15", "20")

expr_file <- "/proj/c.zihao/work1/01data/exprset/LUAD_exprSet_filtered.csv"
expr <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)

sample_info <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
sample_info$group <- substr(sample_info$sample, 14, 16)
tumor_samples <- sample_info$sample[sample_info$group == "01A"]

expr_tumor <- expr[, tumor_samples, drop = FALSE]
gene_sd <- apply(expr_tumor, 1, sd, na.rm = TRUE)

normal_file <- "/proj/c.zihao/work1/01data/03TCGAnormal/paired_normals/TCGA-LUAD_paired_normal.csv"
paired_normal <- read.csv(normal_file, header = TRUE, row.names = 1, check.names = FALSE)
missing_normal_genes <- setdiff(rownames(expr_tumor), rownames(paired_normal))
if (length(missing_normal_genes) > 0L) {
  stop("Paired-normal matrix is missing ", length(missing_normal_genes), " genes from the LUAD tumour expression matrix.")
}
paired_normal <- round(paired_normal[rownames(expr_tumor), , drop = FALSE], 5)

cat("Original samples:", ncol(expr), "\n")
cat("Tumor samples used:", ncol(expr_tumor), "\n")
cat("Paired LUAD normal reference samples:", ncol(paired_normal), "\n")

for (level in LEVELS) {
  noise_level <- as.numeric(level) / 100
  data_dir <- file.path(base_dir, level, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(123)

  noise_mat <- matrix(
    rnorm(
      nrow(expr_tumor) * ncol(expr_tumor),
      mean = 0,
      sd = rep(noise_level * gene_sd, ncol(expr_tumor))
    ),
    nrow = nrow(expr_tumor),
    ncol = ncol(expr_tumor),
    byrow = FALSE
  )
  rownames(noise_mat) <- rownames(expr_tumor)
  colnames(noise_mat) <- colnames(expr_tumor)

  expr_noise <- expr_tumor + noise_mat
  expr_noise[expr_noise < 0] <- 0
  expr_noise <- round(expr_noise, 5)

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

  write.csv(
    paired_normal,
    file = file.path(data_dir, "LUAD_normal_reference.csv"),
    quote = FALSE
  )

  write.csv(
    data.frame(sample = colnames(paired_normal), stringsAsFactors = FALSE),
    file = file.path(data_dir, "LUAD_normal_reference_samples.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  cat(
    "Noise level ", level, "% (sd factor ", noise_level, ")",
    " -> ", data_dir, "\n",
    sep = ""
  )
}
