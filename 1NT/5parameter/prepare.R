rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


base_dir <- paste0(PROJ_ROOT, "/1NT/5parameter")
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- paste0(PROJ_ROOT, "/1NT/1data/exprset/LUAD_exprSet_filtered.csv")
expr <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)

sample_ids <- colnames(expr)
patient_ids <- substr(sample_ids, 1, 12)
sample_types <- substr(sample_ids, 14, 16)

tumor_patients <- unique(patient_ids[sample_types == "01A"])
normal_patients <- unique(patient_ids[sample_types == "11A"])
paired_patients <- sort(intersect(tumor_patients, normal_patients))

paired_samples <- sample_ids[patient_ids %in% paired_patients & sample_types %in% c("01A", "11A")]
paired_samples <- paired_samples[order(substr(paired_samples, 1, 12), substr(paired_samples, 14, 16))]

expr_paired <- expr[, paired_samples, drop = FALSE]

metadata <- data.frame(
  sample = paired_samples,
  patient = substr(paired_samples, 1, 12),
  group = substr(paired_samples, 14, 16),
  stringsAsFactors = FALSE
)
metadata$group_label <- ifelse(metadata$group == "01A", "Tumor", "Normal")

write.csv(expr_paired, file.path(data_dir, "LUAD_paired_expr.csv"), quote = FALSE)
write.csv(metadata, file.path(data_dir, "LUAD_paired_metadata.csv"), row.names = FALSE, quote = FALSE)

cat("Total LUAD samples:", ncol(expr), "\n")
cat("Paired patients:", length(paired_patients), "\n")
cat("Paired samples:", ncol(expr_paired), "\n")
cat("Saved paired expression matrix to:", file.path(data_dir, "LUAD_paired_expr.csv"), "\n")
