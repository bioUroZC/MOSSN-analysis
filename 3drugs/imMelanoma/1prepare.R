# ===================================================

# ===================================================

rm(list = ls())

library(dplyr)

# =======================================================

# =======================================================

dataset_paths <- list(
  GSE78220 = "/proj/c.zihao/work1/3drugs/imMelanoma/GSE78220/data/exprSet.csv",
  GSE91061 = "/proj/c.zihao/work1/3drugs/imMelanoma/GSE91061/data/exprSet.csv",
  GSE100797 = "/proj/c.zihao/work1/3drugs/imMelanoma/GSE100797/data/exprSet.csv",
  Nathanson = "/proj/c.zihao/work1/3drugs/imMelanoma/Nathanson/data/exprSet.csv",
  PRJEB23709 = "/proj/c.zihao/work1/3drugs/imMelanoma/PRJEB23709/data/exprSet.csv"
)


links = read.csv("/proj/c.zihao/work1/1NT/1data/string/links.csv", row.names = 1)
common_genes <- unique(union(links$protein1, links$protein2))
message("Number of common genes: ", length(common_genes))

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
