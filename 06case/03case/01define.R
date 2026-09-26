rm(list = ls())

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(base_dir, "01data/exprSet_filtered.csv")
matrix_file <- file.path(base_dir, "04results/MOSSN/merged_matrix.csv")

module_list <- list(
  IFNG_response = c(
    "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "STAT2",
    "IRF1", "IRF7", "IRF9", "CXCL9", "CXCL10", "CXCL11", "GBP1",
    "GBP2", "GBP5", "IDO1", "ISG15", "IFIT1", "IFIT2", "IFIT3",
    "OAS1", "OAS2", "MX1", "TAP1", "TAP2", "PSMB8", "PSMB9"
  ),
  Antigen_presentation = c(
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M",
    "TAP1", "TAPBP", "ERAP1", "ERAP2", "PSMB8", "PSMB9",
    "PSMB10", "CALR", "CANX", "PDIA3", "NLRC5", "HLA-DRA", "HLA-DRB1",
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA"
  ),
  Cytotoxicity = c(
    "CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1",
    "GNLY", "NKG7", "CTSW", "KLRD1", "KLRK1", "NCR1", "NCR3", "FGFBP2",
    "CCL5", "CX3CR1", "TBX21", "EOMES", "IFNG", "LAMP1", "FASLG", "TNFSF10"
  )
)

checkpoint_genes <- c(
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD80", "CD86", "LAG3", "TIGIT",
  "HAVCR2", "CD28", "ICOS", "ICOSLG", "TNFRSF9", "TNFRSF4", "CD40",
  "CD40LG", "BTLA", "VSIR", "ENTPD1", "LAIR1", "CD27", "CD70"
)

immune_union_genes <- sort(unique(unlist(module_list)))
partner_genes <- setdiff(immune_union_genes, checkpoint_genes)

gene_membership <- data.frame()
for (module_name in names(module_list)) {
  genes <- unique(module_list[[module_name]])
  gene_membership <- rbind(gene_membership, data.frame(
    module = module_name, gene = genes, gene_class = "core_module"
  ))
}
gene_membership <- rbind(
  gene_membership,
  data.frame(module = "Checkpoint_neighborhood", gene = checkpoint_genes,
             gene_class = "checkpoint_seed"),
  data.frame(module = "Checkpoint_neighborhood", gene = partner_genes,
             gene_class = "immune_partner")
)

gene_membership$gene_norm <- gsub("[-._]", "", toupper(gene_membership$gene))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_genes <- rownames(expr)
expr_norm <- gsub("[-._]", "", toupper(expr_genes))

idx <- match(gene_membership$gene_norm, expr_norm)
gene_membership$gene_expr <- expr_genes[idx]

edge_names <- data.table::fread(matrix_file, select = "Interaction")[[1]]
edge_names <- unique(edge_names)
edge_names <- edge_names[grepl("_", edge_names, fixed = TRUE)]

edge_df <- data.frame(
  edge = edge_names,
  gene1 = sub("_.*$", "", edge_names),
  gene2 = sub("^[^_]*_", "", edge_names)
)
edge_df <- edge_df[edge_df$gene1 != "" & edge_df$gene2 != "", ]
edge_df$gene1_norm <- gsub("[-._]", "", toupper(edge_df$gene1))
edge_df$gene2_norm <- gsub("[-._]", "", toupper(edge_df$gene2))

backbone_genes <- c(edge_df$gene1, edge_df$gene2)
backbone_norm <- c(edge_df$gene1_norm, edge_df$gene2_norm)

idx <- match(gene_membership$gene_norm, backbone_norm)
gene_membership$gene_backbone <- backbone_genes[idx]

gene_membership$in_expression <- !is.na(gene_membership$gene_expr)
gene_membership$in_backbone <- !is.na(gene_membership$gene_backbone)

module_edges <- data.frame()
for (module_name in names(module_list)) {
  genes <- gsub("[-._]", "", toupper(module_list[[module_name]]))
  keep <- edge_df$gene1_norm %in% genes & edge_df$gene2_norm %in% genes
  one <- edge_df[keep, c("edge", "gene1", "gene2")]
  one$module <- rep(module_name, nrow(one))
  one$edge_type <- rep("within_module", nrow(one))
  module_edges <- rbind(module_edges, one)
}

checkpoint_norm <- gsub("[-._]", "", toupper(checkpoint_genes))
immune_norm <- gsub("[-._]", "", toupper(immune_union_genes))
keep <- (edge_df$gene1_norm %in% checkpoint_norm & edge_df$gene2_norm %in% immune_norm) |
  (edge_df$gene2_norm %in% checkpoint_norm & edge_df$gene1_norm %in% immune_norm)
one <- edge_df[keep, c("edge", "gene1", "gene2")]
one$module <- rep("Checkpoint_neighborhood", nrow(one))
one$edge_type <- rep("checkpoint_neighborhood", nrow(one))
module_edges <- rbind(module_edges, one)

module_edges <- module_edges[, c("module", "edge", "gene1", "gene2", "edge_type")]
module_edges <- unique(module_edges)

write.csv(gene_membership, file.path(out_dir, "module_gene_membership.csv"),
          row.names = FALSE)
write.csv(module_edges, file.path(out_dir, "module_edge_membership.csv"),
          row.names = FALSE)
write.csv(gene_membership[, c("module", "gene")],
          file.path(out_dir, "module_gene_list.csv"),
          row.names = FALSE, quote = FALSE)

module_order <- c("IFNG_response", "Antigen_presentation",
                  "Cytotoxicity", "Checkpoint_neighborhood")
module_qc <- data.frame()
for (module_name in module_order) {
  g <- gene_membership[gene_membership$module == module_name, ]
  e <- module_edges[module_edges$module == module_name, ]
  module_qc <- rbind(module_qc, data.frame(
    module = module_name,
    defined_genes = length(unique(g$gene)),
    mapped_expression_genes = length(unique(g$gene[g$in_expression])),
    mapped_backbone_genes = length(unique(g$gene[g$in_backbone])),
    mapped_edges = nrow(e),
    mapped_genes_in_edges = length(unique(c(e$gene1, e$gene2)))
  ))
}
print(module_qc)
