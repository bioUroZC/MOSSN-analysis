rm(list = ls())
dataset_name <- "HNSC"
source("/proj/c.zihao/work1/function/benchmark/LIONESS.R")

base_dir <- "/proj/c.zihao/work1/02analysis/03net/2HuRI/out"
save_path <- file.path(base_dir, "LIONESS", dataset_name)

expr_file <- file.path(
  "/proj/c.zihao/work1/01data/exprset",
  paste0(dataset_name, "_exprSet_filtered.csv")
)
ppi_file <- "/proj/c.zihao/work1/01data/05HuRI/huri_link.csv"

unlink(save_path, recursive = TRUE, force = TRUE)
dir.create(save_path, recursive = TRUE)
write.csv(
  LIONcal(expr_file, ppi_file),
  file.path(save_path, "result.csv"),
  row.names = FALSE
)
