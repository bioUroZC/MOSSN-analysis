args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "/proj/c.zihao/work1/01data/03TCGA"
output_dir <- if (length(args) >= 2) args[[2]] else "/proj/c.zihao/work1/01data/03TCGAnormal"
normal_dir <- file.path(output_dir, "paired_normals")

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}
dir.create(normal_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(
  input_dir,
  pattern = "^TCGA-[A-Z]+\\.csv$",
  full.names = TRUE
)
if (length(input_files) == 0L) {
  stop("No TCGA expression matrices matching '^TCGA-[A-Z]+\\.csv$' found in ", input_dir)
}

parse_barcode <- function(sample_ids) {
  fields <- strsplit(sample_ids, "_", fixed = TRUE)
  valid <- lengths(fields) >= 4L
  patient <- vapply(fields, function(x) {
    if (length(x) < 4L) return(NA_character_)
    paste(x[1:3], collapse = "_")
  }, character(1))
  aliquot <- vapply(fields, function(x) {
    if (length(x) < 4L) return(NA_character_)
    x[[4]]
  }, character(1))

  data.frame(
    sample = sample_ids,
    patient = patient,
    sample_type = substr(aliquot, 1L, 2L),
    valid = valid,
    stringsAsFactors = FALSE
  )
}

manifests <- list()
summaries <- list()

for (input_file in sort(input_files)) {
  cancer <- sub("^TCGA-([A-Z]+)\\.csv$", "\\1", basename(input_file))
  message("Processing TCGA-", cancer)

  expression <- read.csv(
    input_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_info <- parse_barcode(colnames(expression))

  if (any(!sample_info$valid)) {
    warning(
      "Ignoring malformed sample barcodes in TCGA-", cancer, ": ",
      paste(sample_info$sample[!sample_info$valid], collapse = ", ")
    )
    sample_info <- sample_info[sample_info$valid, , drop = FALSE]
  }

  tumour <- sample_info[sample_info$sample_type == "01", c("patient", "sample"), drop = FALSE]
  normal <- sample_info[sample_info$sample_type == "11", c("patient", "sample"), drop = FALSE]
  tumour <- tumour[order(tumour$patient, tumour$sample), , drop = FALSE]
  normal <- normal[order(normal$patient, normal$sample), , drop = FALSE]
  tumour <- tumour[!duplicated(tumour$patient), , drop = FALSE]
  normal <- normal[!duplicated(normal$patient), , drop = FALSE]
  names(tumour)[2] <- "tumour_sample"
  names(normal)[2] <- "normal_sample"

  pairs <- merge(tumour, normal, by = "patient", all = FALSE, sort = TRUE)
  pairs$cancer <- cancer
  pairs <- pairs[, c("cancer", "patient", "tumour_sample", "normal_sample")]

  normal_matrix <- expression[, pairs$normal_sample, drop = FALSE]
  output_file <- file.path(normal_dir, paste0("TCGA-", cancer, "_paired_normal.csv"))
  write.csv(normal_matrix, output_file, quote = FALSE)

  manifests[[cancer]] <- pairs
  summaries[[cancer]] <- data.frame(
    cancer = cancer,
    input_samples = ncol(expression),
    primary_tumour_samples = sum(sample_info$sample_type == "01"),
    solid_tissue_normal_samples = sum(sample_info$sample_type == "11"),
    paired_patients = nrow(pairs),
    output_file = output_file,
    stringsAsFactors = FALSE
  )

  message(
    "  paired patients: ", nrow(pairs),
    "; normal matrix: ", nrow(normal_matrix), " genes x ", ncol(normal_matrix), " samples"
  )
}

manifest <- do.call(rbind, manifests)
summary_table <- do.call(rbind, summaries)
write.csv(manifest, file.path(output_dir, "paired_sample_manifest.csv"), row.names = FALSE, quote = FALSE)
write.csv(summary_table, file.path(output_dir, "extraction_summary.csv"), row.names = FALSE, quote = FALSE)

message("\nWrote ", nrow(manifest), " paired normal samples from ", nrow(summary_table), " TCGA projects.")
