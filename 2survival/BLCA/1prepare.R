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

# =======================================================

dataset_paths <- list(
  GSE13507 = paste0(PROJ_ROOT, "/2survival/BLCA/GSE13507/data/exprSet.csv"),
  GSE31684 = paste0(PROJ_ROOT, "/2survival/BLCA/GSE31684/data/exprSet.csv"),
  GSE32894 = paste0(PROJ_ROOT, "/2survival/BLCA/GSE32894/data/exprSet.csv"),
  GSE48276 = paste0(PROJ_ROOT, "/2survival/BLCA/GSE48276/data/exprSet.csv"),
  TCGABLCA = paste0(PROJ_ROOT, "/2survival/BLCA/TCGABLCA/data/exprSet.csv")
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
