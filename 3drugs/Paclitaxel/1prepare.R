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
  GSE194040 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE194040/data/exprSet.csv",
  GSE20194 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE20194/data/exprSet.csv",
  GSE20271 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE20271/data/exprSet.csv",
  
  GSE241876 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE241876/data/exprSet.csv",
  GSE28844 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE28844/data/exprSet.csv",
  GSE32646 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE32646/data/exprSet.csv",
  GSE41998 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE41998/data/exprSet.csv",
  
  GSE50948 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE50948/data/exprSet.csv",
  GSE66305 = "/proj/c.zihao/work1/3drugs/Paclitaxel/GSE66305/data/exprSet.csv"
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
  
  print(min(expr_filtered))
  print(max(expr_filtered))
  
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
