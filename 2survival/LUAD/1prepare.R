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
  GSE30219 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE30219/data/exprSet.csv"),

  GSE31210 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE31210/data/exprSet.csv"),
  GSE41271 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE41271/data/exprSet.csv"),
  GSE42127 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE42127/data/exprSet.csv"),
  GSE50081 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE50081/data/exprSet.csv"),
  
  GSE68465 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE68465/data/exprSet.csv"),
  GSE72094 = paste0(PROJ_ROOT, "/2survival/LUAD/GSE72094/data/exprSet.csv"),
  TCGALUAD = paste0(PROJ_ROOT, "/2survival/LUAD/TCGALUAD/data/exprSet.csv")
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
