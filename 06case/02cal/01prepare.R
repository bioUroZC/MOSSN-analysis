rm(list = ls())

expr_file <- "/proj/c.zihao/work1/06case/01data/exprSet.csv"
links_file <- "/proj/c.zihao/work1/01data/06string/links.csv"
out_file <- "/proj/c.zihao/work1/06case/01data/exprSet_filtered.csv"

links <- read.csv(links_file, row.names = 1)
common_genes <- unique(union(links$protein1, links$protein2))
message("STRING genes: ", length(common_genes))

expr <- read.csv(expr_file, row.names = 1)
expr_filtered <- expr[rownames(expr) %in% common_genes, ]
expr_filtered <- log2(expr_filtered + 1)

write.csv(expr_filtered, file = out_file)
message("Saved ", nrow(expr_filtered), " genes x ", ncol(expr_filtered),
        " samples to ", out_file)
