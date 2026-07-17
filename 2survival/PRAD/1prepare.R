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
  DKFZ2018 = paste0(PROJ_ROOT, "/2survival/PRAD/DKFZ2018/data/exprSet.csv"),
  GSE116918 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE116918/data/exprSet.csv"),
  GSE21034 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE21034/data/exprSet.csv"),

  GSE46602 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE46602/data/exprSet.csv"),
  GSE54460 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE54460/data/exprSet.csv"),
  GSE70768 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE70768/data/exprSet.csv"),
  GSE70769 = paste0(PROJ_ROOT, "/2survival/PRAD/GSE70769/data/exprSet.csv"),

  TCGAPRAD = paste0(PROJ_ROOT, "/2survival/PRAD/TCGAPRAD/data/exprSet.csv")
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
