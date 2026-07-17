# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

# =======================================================
# Clear Environment & Load Libraries
# =======================================================

dataset_paths <- list(
  ACICAM = paste0(PROJ_ROOT, "/2survival/CRC/ACICAM/data/exprSet.csv"),
  GSE12945 = paste0(PROJ_ROOT, "/2survival/CRC/GSE12945/data/exprSet.csv"),
  GSE17536 = paste0(PROJ_ROOT, "/2survival/CRC/GSE17536/data/exprSet.csv"),
  GSE17537 = paste0(PROJ_ROOT, "/2survival/CRC/GSE17537/data/exprSet.csv"),
  GSE28722 = paste0(PROJ_ROOT, "/2survival/CRC/GSE28722/data/exprSet.csv"),
  
  GSE29621 = paste0(PROJ_ROOT, "/2survival/CRC/GSE29621/data/exprSet.csv"),
  GSE39582 = paste0(PROJ_ROOT, "/2survival/CRC/GSE39582/data/exprSet.csv"),
  GSE41258 = paste0(PROJ_ROOT, "/2survival/CRC/GSE41258/data/exprSet.csv"),
  TCGACRC = paste0(PROJ_ROOT, "/2survival/CRC/TCGACRC/data/exprSet.csv")
  
)

links = read.csv(paste0(PROJ_ROOT, "/1NT/1data/string/links.csv"), row.names = 1)
common_genes <- unique(union(links$protein1, links$protein2))
message("Number of common genes: ", length(common_genes))

names(dataset_paths)

# =======================================================

# =======================================================


for (name in names(dataset_paths)) {
  
  original_path <- dataset_paths[[name]]
  expr <- read.csv(original_path, row.names = 1)
  
  expr_filtered <- expr[which(rownames(expr) %in% common_genes), ]
  
  expr_filtered <- log2(expr_filtered + 1) 
  
  print(min(expr_filtered))
  print(max(expr_filtered))
  
  dir_path <- dirname(original_path)
  base_name <- tools::file_path_sans_ext(basename(original_path))
  save_path <- file.path(dir_path, paste0(base_name, "_filtered.csv"))
  
  write.csv(expr_filtered, file = save_path)
  
  print(dim(expr_filtered))
  
  message("Saved filtered file to: ", save_path)
}
