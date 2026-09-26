rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

dataset_paths <- list(
GSE11121 = "/proj/c.zihao/work1/03survival/BRCA/GSE11121/data/exprSet.csv",
GSE162228 = "/proj/c.zihao/work1/03survival/BRCA/GSE162228/data/exprSet.csv",
GSE17705 = "/proj/c.zihao/work1/03survival/BRCA/GSE17705/data/exprSet.csv",
GSE20685 = "/proj/c.zihao/work1/03survival/BRCA/GSE20685/data/exprSet.csv",

GSE20711 = "/proj/c.zihao/work1/03survival/BRCA/GSE20711/data/exprSet.csv",
GSE21653 = "/proj/c.zihao/work1/03survival/BRCA/GSE21653/data/exprSet.csv",
GSE22219 = "/proj/c.zihao/work1/03survival/BRCA/GSE22219/data/exprSet.csv",
GSE25055 = "/proj/c.zihao/work1/03survival/BRCA/GSE25055/data/exprSet.csv",
GSE25065 = "/proj/c.zihao/work1/03survival/BRCA/GSE25065/data/exprSet.csv",

GSE45255 = "/proj/c.zihao/work1/03survival/BRCA/GSE45255/data/exprSet.csv",
GSE61304 = "/proj/c.zihao/work1/03survival/BRCA/GSE61304/data/exprSet.csv",
GSE7390 = "/proj/c.zihao/work1/03survival/BRCA/GSE7390/data/exprSet.csv",

TCGABRCA = "/proj/c.zihao/work1/03survival/BRCA/TCGABRCA/data/exprSet.csv"

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
