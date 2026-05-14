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
  GSE13861 = "/proj/c.zihao/work1/2survival/STAD/GSE13861/data/exprSet.csv",
  GSE15459 = "/proj/c.zihao/work1/2survival/STAD/GSE15459/data/exprSet.csv",
  GSE26253 = "/proj/c.zihao/work1/2survival/STAD/GSE26253/data/exprSet.csv",
  GSE26899 = "/proj/c.zihao/work1/2survival/STAD/GSE26899/data/exprSet.csv",
  GSE26901 = "/proj/c.zihao/work1/2survival/STAD/GSE26901/data/exprSet.csv",
  GSE29272 = "/proj/c.zihao/work1/2survival/STAD/GSE29272/data/exprSet.csv",
  GSE57303 = "/proj/c.zihao/work1/2survival/STAD/GSE57303/data/exprSet.csv",
  GSE62254 = "/proj/c.zihao/work1/2survival/STAD/GSE62254/data/exprSet.csv",
  GSE84437 = "/proj/c.zihao/work1/2survival/STAD/GSE84437/data/exprSet.csv",
  TCGASTAD = "/proj/c.zihao/work1/2survival/STAD/TCGASTAD/data/exprSet.csv"
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
