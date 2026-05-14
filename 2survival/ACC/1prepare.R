# ===================================================

# ===================================================

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

# =======================================================

# =======================================================

dataset_paths <- list(
  GSE10927 = "/proj/c.zihao/work1/2survival/ACC/GSE10927/data/exprSet.csv",
  GSE19750 = "/proj/c.zihao/work1/2survival/ACC/GSE19750/data/exprSet.csv",
  GSE33371 = "/proj/c.zihao/work1/2survival/ACC/GSE33371/data/exprSet.csv",
  GSE49278 = "/proj/c.zihao/work1/2survival/ACC/GSE49278/data/exprSet.csv",
  GSE76019 = "/proj/c.zihao/work1/2survival/ACC/GSE76019/data/exprSet.csv",
  GSE76021 = "/proj/c.zihao/work1/2survival/ACC/GSE76021/data/exprSet.csv",
  TCGAACC = "/proj/c.zihao/work1/2survival/ACC/TCGAACC/data/exprSet.csv"
)


links = read.csv("/proj/c.zihao/work1/1NT/1data/string/links.csv", row.names = 1)
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
