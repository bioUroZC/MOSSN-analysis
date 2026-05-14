rm(list=ls())

cancertypes <- c("BLCA", "LIHC", "LUAD", "SARC", "STAD")

source_root <- "/proj/c.zihao/work1/2survival"
output_dir <- "/proj/c.zihao/work1/4coupled/EXPout"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (tumor in cancertypes) {
  source_file <- file.path(
    source_root,
    tumor,
    paste0("TCGA", tumor),
    "data",
    "exprSet_filtered.csv"
  )

  if (!file.exists(source_file)) {
    warning(sprintf("File not found: %s", source_file))
    next
  }

  target_file <- file.path(
    output_dir,
    sprintf("EXP_%s.csv", tumor)
  )

  file.copy(source_file, target_file, overwrite = TRUE)
  message(sprintf("Copied %s -> %s", source_file, target_file))
}
