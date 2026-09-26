rm(list = ls())

base_dir <- "/proj/c.zihao/work1/02analysis/01robust"
LEVELS <- c(10, 30, 50, 70)

expr_file <- "/proj/c.zihao/work1/01data/exprset/LUAD_exprSet_filtered.csv"
expr <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)

sample_info <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
sample_info$group <- substr(sample_info$sample, 14, 16)
tumor_samples <- sample_info$sample[sample_info$group == "01A"]

normal_file <- "/proj/c.zihao/work1/01data/03TCGAnormal/paired_normals/TCGA-LUAD_paired_normal.csv"
paired_normal <- read.csv(normal_file, header = TRUE, row.names = 1, check.names = FALSE)
missing_normal_genes <- setdiff(rownames(expr), rownames(paired_normal))
if (length(missing_normal_genes) > 0L) {
  stop("Paired-normal matrix is missing ", length(missing_normal_genes), " genes from the LUAD expression matrix.")
}
paired_normal <- paired_normal[rownames(expr), , drop = FALSE]
normal_n <- ncol(paired_normal)

cat("Original samples:", ncol(expr), "\n")
cat("Original tumor samples:", length(tumor_samples), "\n")
cat("Original paired LUAD normal samples:", normal_n, "\n")

for (level in LEVELS) {
  ratio <- level / 100
  data_dir <- file.path(base_dir, level, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(123)

  selected_samples <- sort(
    sample(
      tumor_samples,
      size = max(1, ceiling(length(tumor_samples) * ratio)),
      replace = FALSE
    )
  )
  expr_half <- expr[, selected_samples, drop = FALSE]
  expr_half <- round(expr_half, 5)

  write.csv(
    expr_half,
    file = file.path(data_dir, "LUAD_exprSet_half.csv"),
    quote = FALSE
  )

  write.csv(
    data.frame(
      sample = selected_samples,
      group = substr(selected_samples, 14, 16),
      stringsAsFactors = FALSE
    ),
    file = file.path(data_dir, "LUAD_half_samples.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  normal_half_n <- max(3, ceiling(normal_n * ratio))
  normal_keep <- sort(sample(colnames(paired_normal), size = normal_half_n, replace = FALSE))
  normal_half <- round(paired_normal[, normal_keep, drop = FALSE], 5)

  write.csv(
    normal_half,
    file = file.path(data_dir, "LUAD_normal_half.csv"),
    quote = FALSE
  )

  write.csv(
    data.frame(sample = normal_keep, stringsAsFactors = FALSE),
    file = file.path(data_dir, "LUAD_normal_half_samples.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  cat(
    "Selected ", level, "% tumor samples: ", length(selected_samples),
    ", paired LUAD normal samples: ", normal_half_n,
    " -> ", data_dir, "\n",
    sep = ""
  )
}
