# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(data.table)
library(tidyr)

# =======================================================

# =======================================================


dataset_paths <- list(
  BLCA = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-BLCA.csv"),
  BRCA = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-BRCA.csv"),
  CRC = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-CRC.csv"),
  ESCA = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-ESCA.csv"),
  
  HNSC = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-HNSC.csv"),
  KIRC = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-KIRC.csv"),
  LIHC = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-LIHC.csv"),
  LUAD = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-LUAD.csv"),
  
  LUSC = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-LUSC.csv"),
  PRAD = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-PRAD.csv"),
  STAD = paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-STAD.csv")
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
  print(expr_filtered[1:4,1:4])
  
  expr_filtered <- log2(expr_filtered + 1) 
  
  print(min(expr_filtered))
  print(max(expr_filtered))
  
  dir_path <- dirname(original_path)
  save_path <- paste0(paste0(PROJ_ROOT, '/1NT/1data/exprset/'), name, "_exprSet_filtered.csv")
  
  write.csv(expr_filtered, file = save_path)
  
  print(dim(expr_filtered))
  
  message("Saved filtered file to: ", save_path)
}
