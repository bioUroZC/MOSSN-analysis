rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

dataset_paths <- list(
  GSE102073 = "/proj/c.zihao/work1/03survival/OV/GSE102073/data/exprSet.csv",
  GSE26193  = "/proj/c.zihao/work1/03survival/OV/GSE26193/data/exprSet.csv",
  GSE26712  = "/proj/c.zihao/work1/03survival/OV/GSE26712/data/exprSet.csv",
  GSE30161  = "/proj/c.zihao/work1/03survival/OV/GSE30161/data/exprSet.csv",
  GSE31245  = "/proj/c.zihao/work1/03survival/OV/GSE31245/data/exprSet.csv",
  GSE32062  = "/proj/c.zihao/work1/03survival/OV/GSE32062/data/exprSet.csv",
  GSE51088  = "/proj/c.zihao/work1/03survival/OV/GSE51088/data/exprSet.csv",
  GSE53963  = "/proj/c.zihao/work1/03survival/OV/GSE53963/data/exprSet.csv",
  GSE8842   = "/proj/c.zihao/work1/03survival/OV/GSE8842/data/exprSet.csv",
  GSE9891   = "/proj/c.zihao/work1/03survival/OV/GSE9891/data/exprSet.csv",
  MTAB386   = "/proj/c.zihao/work1/03survival/OV/MTAB386/data/exprSet.csv",
  TCGAOV    = "/proj/c.zihao/work1/03survival/OV/TCGAOV/data/exprSet.csv"
)

links = read.csv("/proj/c.zihao/work1/01data/06string/links.csv", row.names = 1)
common_genes <- unique(union(links$protein1, links$protein2))
message("Number of common genes: ", length(common_genes))

names(dataset_paths)

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
  
  expr_filtered <- round(expr_filtered, 5)
  
  write.csv(expr_filtered, file = save_path)
  
  print(dim(expr_filtered))
  
  message("Saved filtered file to: ", save_path)
}
