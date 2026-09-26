rm(list = ls())
library(dplyr)
set.seed(42)

dataset_name <- "LUAD"
base_dir <- "/proj/c.zihao/work1/02analysis/01robust"
LEVELS <- c(10, 30, 50, 70)
ppiFile <- "/proj/c.zihao/work1/01data/06string/links.csv"

source('/proj/c.zihao/work1/function/benchmark/LIONESS.R')

for (level in LEVELS) {
  cat("#==========", dataset_name, paste0(level, "%"), "==========\n")
  level_dir <- file.path(base_dir, level)
  save_path <- file.path(level_dir, "LIONESS", dataset_name)
  exprSetFile <- file.path(level_dir, "data", "LUAD_exprSet_half.csv")

  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

  resultDF <- LIONcal(exprSetFile, ppiFile)
  num_cols <- vapply(resultDF, is.numeric, logical(1))
  resultDF[num_cols] <- round(resultDF[num_cols], 5)

  write.csv(resultDF, file = file.path(save_path, "result.csv"), row.names = FALSE)
}
