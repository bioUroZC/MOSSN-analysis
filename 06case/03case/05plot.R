rm(list = ls())

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")

expr_file <- file.path(base_dir, "01data/exprSet_filtered.csv")
matrix_file <- file.path(base_dir, "04results/MOSSN/merged_matrix.csv")

case_sample <- "SAM2eb07dedf07f"

checkpoint_genes <- c(
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD80", "CD86", "LAG3", "TIGIT",
  "HAVCR2", "CD28", "ICOS", "ICOSLG", "TNFRSF9", "TNFRSF4", "CD40",
  "CD40LG", "BTLA", "VSIR", "ENTPD1", "LAIR1", "CD27", "CD70"
)

keep_edges <- c(
  "B2M_CD8A",
  "CD274_CD8A", "CD8A_PDCD1", "CD274_PDCD1",
  "CD86_CTLA4", "CD80_CTLA4", "CD274_CTLA4",
  "CD28_CD8A", "CD28_CD86", "CD28_CD80",
  "CD8A_FASLG", "CD8A_IFNG",
  "CD8A_LAG3",
  "CD27_CD8A", "CD27_CD40",
  "CD40_CD86", "CD40_CD80"
)

link <- data.table::fread(matrix_file, select = c("Interaction", case_sample),
                          data.table = FALSE)
edge_df <- data.frame(
  gene1 = sub("_.*$", "", link$Interaction),
  gene2 = sub("^[^_]*_", "", link$Interaction),
  weight = link[[case_sample]],
  edge = link$Interaction
)
edge_df <- edge_df[edge_df$edge %in% keep_edges & edge_df$weight > 0, ]
edge_df <- edge_df[order(edge_df$weight, decreasing = TRUE), ]
edge_df$edge_type <- "checkpoint_core"

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)

node_df <- data.frame(gene = sort(unique(c(edge_df$gene1, edge_df$gene2))))
node_expr <- as.matrix(expr[node_df$gene, ])
node_z <- (node_expr - rowMeans(node_expr)) / apply(node_expr, 1, sd)
node_df$expr_z <- node_z[, case_sample]
node_df$role <- ifelse(node_df$gene %in% checkpoint_genes, "Checkpoint gene", "Immune partner")

idx1 <- match(edge_df$gene1, node_df$gene)
idx2 <- match(edge_df$gene2, node_df$gene)
cytoscape_df <- data.frame(
  gene1 = edge_df$gene1,
  gene2 = edge_df$gene2,
  weight = round(edge_df$weight, 3),
  edge_type = edge_df$edge_type,
  gene1_expr_z = round(node_df$expr_z[idx1], 3),
  gene1_role = node_df$role[idx1],
  gene2_expr_z = round(node_df$expr_z[idx2], 3),
  gene2_role = node_df$role[idx2]
)

node_df$expr_z <- round(node_df$expr_z, 5)

write.csv(cytoscape_df, file.path(out_dir, "cytoscape_network.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(node_df, file.path(out_dir, "cytoscape_node_attr.csv"),
          row.names = FALSE, quote = FALSE)
