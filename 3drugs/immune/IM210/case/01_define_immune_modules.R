rm(list = ls())

library(data.table)
library(dplyr)

base_dir <- "/proj/c.zihao/work1/3drugs"
code_dir <- file.path(base_dir, "immune/IM210/case")
out_dir <- file.path(base_dir, "immune/IM210/case/results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_gene_symbol <- function(x) {
  x <- toupper(x)
  gsub("[-._]", "", x)
}

# Predefined immune modules are fixed before any response analysis.
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

gene_membership <- bind_rows(
  lapply(names(module_list), function(module_name) {
    data.frame(
      module = module_name,
      gene = unique(module_list[[module_name]]),
      gene_class = "core_module",
      stringsAsFactors = FALSE
    )
  }),
  data.frame(
    module = "Checkpoint_neighborhood",
    gene = unique(checkpoint_genes),
    gene_class = "checkpoint_seed",
    stringsAsFactors = FALSE
  ),
  data.frame(
    module = "Checkpoint_neighborhood",
    gene = setdiff(immune_union_genes, checkpoint_genes),
    gene_class = "immune_partner",
    stringsAsFactors = FALSE
  )
) %>%
  distinct(module, gene, gene_class) %>%
  mutate(gene_norm = normalize_gene_symbol(gene))

cat("Step 1: load expression genes\n")
expr_path <- file.path(base_dir, "immune/IM210/data/exprSet_filtered.csv")
expr_genes <- fread(expr_path, select = 1, data.table = FALSE)[[1]]
expr_genes <- unique(expr_genes)
expr_gene_map <- data.frame(
  gene_expr = expr_genes,
  gene_norm = normalize_gene_symbol(expr_genes),
  stringsAsFactors = FALSE
) %>%
  distinct(gene_norm, .keep_all = TRUE)

cat("Step 2: load backbone links\n")
link_path <- file.path(base_dir, "immune/1Matrix/MOSSN_noPrior_Matrix.csv")
link_dt <- fread(link_path, select = "Interaction", data.table = FALSE)
link_vec <- unique(link_dt$Interaction)

edge_df <- data.frame(
  edge = link_vec,
  stringsAsFactors = FALSE
) %>%
  filter(grepl("_", edge, fixed = TRUE)) %>%
  tidyr::separate(
    col = edge,
    into = c("gene1", "gene2"),
    sep = "_",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  filter(!is.na(gene1), !is.na(gene2), gene1 != "", gene2 != "") %>%
  mutate(
    gene1_norm = normalize_gene_symbol(gene1),
    gene2_norm = normalize_gene_symbol(gene2)
  )

backbone_gene_map <- bind_rows(
  data.frame(gene_backbone = edge_df$gene1, gene_norm = edge_df$gene1_norm, stringsAsFactors = FALSE),
  data.frame(gene_backbone = edge_df$gene2, gene_norm = edge_df$gene2_norm, stringsAsFactors = FALSE)
) %>%
  distinct(gene_norm, .keep_all = TRUE)

cat("Step 3: map genes and edges to predefined modules\n")
within_module_edges <- bind_rows(
  lapply(names(module_list), function(module_name) {
    genes <- unique(normalize_gene_symbol(module_list[[module_name]]))
    edge_df %>%
      filter(gene1_norm %in% genes, gene2_norm %in% genes) %>%
      mutate(
        module = module_name,
        edge_type = "within_module"
      )
  })
)

checkpoint_edges <- edge_df %>%
  filter(
    (gene1_norm %in% normalize_gene_symbol(checkpoint_genes) &
       gene2_norm %in% normalize_gene_symbol(immune_union_genes)) |
      (gene2_norm %in% normalize_gene_symbol(checkpoint_genes) &
         gene1_norm %in% normalize_gene_symbol(immune_union_genes))
  ) %>%
  mutate(
    module = "Checkpoint_neighborhood",
    edge_type = "checkpoint_neighborhood"
  )

module_edges <- bind_rows(within_module_edges, checkpoint_edges) %>%
  distinct(module, edge, gene1, gene2, edge_type)

gene_membership <- gene_membership %>%
  left_join(expr_gene_map, by = "gene_norm") %>%
  left_join(backbone_gene_map, by = "gene_norm") %>%
  mutate(
    in_expression = !is.na(gene_expr),
    in_backbone = !is.na(gene_backbone)
  )

edge_qc <- module_edges %>%
  group_by(module) %>%
  summarise(
    mapped_edges = n(),
    mapped_genes_in_edges = n_distinct(c(gene1, gene2)),
    .groups = "drop"
  )

gene_qc <- gene_membership %>%
  group_by(module) %>%
  summarise(
    defined_genes = n_distinct(gene),
    mapped_expression_genes = n_distinct(gene[in_expression]),
    mapped_backbone_genes = n_distinct(gene[in_backbone]),
    .groups = "drop"
  )

module_qc <- gene_qc %>%
  left_join(edge_qc, by = "module") %>%
  mutate(
    mapped_edges = ifelse(is.na(mapped_edges), 0L, mapped_edges),
    mapped_genes_in_edges = ifelse(is.na(mapped_genes_in_edges), 0L, mapped_genes_in_edges)
  ) %>%
  arrange(match(
    module,
    c("IFNG_response", "Antigen_presentation", "Cytotoxicity", "Checkpoint_neighborhood")
  ))

write.csv(
  gene_membership,
  file.path(out_dir, "module_gene_membership.csv"),
  row.names = FALSE
)
write.csv(
  module_edges,
  file.path(out_dir, "module_edge_membership.csv"),
  row.names = FALSE
)
write.csv(
  gene_membership[, c("module", "gene")],
  file.path(out_dir, "module_gene_list.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Finished module definition step.\n")
print(module_qc)
