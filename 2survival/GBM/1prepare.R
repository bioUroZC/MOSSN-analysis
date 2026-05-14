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
  CGGA301 = "/proj/c.zihao/work1/2survival/GBM/CGGA301/data/exprSet.csv",
  CGGA325 = "/proj/c.zihao/work1/2survival/GBM/CGGA325/data/exprSet.csv",
  CGGA693 = "/proj/c.zihao/work1/2survival/GBM/CGGA693/data/exprSet.csv",
  GSE13041 = "/proj/c.zihao/work1/2survival/GBM/GSE13041/data/exprSet.csv",
  GSE16011 = "/proj/c.zihao/work1/2survival/GBM/GSE16011/data/exprSet.csv",
  GSE4412 = "/proj/c.zihao/work1/2survival/GBM/GSE4412/data/exprSet.csv",
  GSE72951 = "/proj/c.zihao/work1/2survival/GBM/GSE72951/data/exprSet.csv",
  GSE74187 = "/proj/c.zihao/work1/2survival/GBM/GSE74187/data/exprSet.csv",
  GSE83300 = "/proj/c.zihao/work1/2survival/GBM/GSE83300/data/exprSet.csv",
  TCGAGBM = "/proj/c.zihao/work1/2survival/GBM/TCGAGBM/data/exprSet.csv"
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
