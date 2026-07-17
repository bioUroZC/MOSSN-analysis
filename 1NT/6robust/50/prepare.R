rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


set.seed(123)

base_dir <- paste0(PROJ_ROOT, "/1NT/6robust/50")
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- paste0(PROJ_ROOT, "/1NT/1data/exprset/LUAD_exprSet_filtered.csv")
expr <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)

sample_info <- data.frame(sample = colnames(expr), stringsAsFactors = FALSE)
sample_info$group <- substr(sample_info$sample, 14, 16)
tumor_samples <- sample_info$sample[sample_info$group == "01A"]
selected_samples <- sort(
  sample(
    tumor_samples,
    size = max(1, ceiling(length(tumor_samples) * 0.5)),
    replace = FALSE
  )
)

sample_num <- ncol(expr)
tumor_num <- length(tumor_samples)
half_num <- length(selected_samples)
expr_half <- expr[, selected_samples, drop = FALSE]

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

normal_file <- paste0(PROJ_ROOT, "/0ref/Test8/combined_expr_df.csv")
normal_df <- read.csv(normal_file, header = TRUE, check.names = FALSE)
lung_normal_df <- normal_df[normal_df$organ == "Lung", , drop = FALSE]

lung_n <- nrow(lung_normal_df)
lung_half_n <- max(1, ceiling(lung_n * 0.5))
lung_keep_idx <- sort(sample(seq_len(lung_n), size = lung_half_n, replace = FALSE))
lung_normal_half <- lung_normal_df[lung_keep_idx, , drop = FALSE]

write.csv(
  lung_normal_half,
  file = file.path(data_dir, "Lung_normal_half.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    sample = lung_normal_half[[1]],
    organ = lung_normal_half$organ,
    stringsAsFactors = FALSE
  ),
  file = file.path(data_dir, "Lung_normal_half_samples.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Original samples:", sample_num, "\n")
cat("Original tumor samples:", tumor_num, "\n")
cat("Selected 50% tumor samples:", half_num, "\n")
cat("Saved tumor-only subset expression matrix to:", file.path(data_dir, "LUAD_exprSet_half.csv"), "\n")
cat("Original Lung normal samples:", lung_n, "\n")
cat("Selected 50% Lung normal samples:", lung_half_n, "\n")
cat("Saved normal reference subset to:", file.path(data_dir, "Lung_normal_half.csv"), "\n")
