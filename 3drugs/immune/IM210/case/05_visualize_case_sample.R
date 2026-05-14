rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(ggraph)

base_dir <- "/proj/c.zihao/work1/3drugs"
out_dir <- file.path(base_dir, "immune/IM210/case/results")

case_sample <- "SAM2eb07dedf07f"
checkpoint_core_genes <- c(
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
  "CD80", "CD86", "CD28", "ICOS", "CD27", "CD40", "CD40LG",
  "CD8A", "IFNG", "B2M", "FASLG"
)

normalize_gene_symbol <- function(x) {
  x <- toupper(x)
  gsub("[-._]", "", x)
}

extract_edge_matrix_subset <- function(matrix_path, edge_ids) {
  edge_ids <- unique(as.character(edge_ids))
  hit_lines <- character(0)

  con <- file(matrix_path, open = "r")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1)

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    parts <- strsplit(line, ",", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      next
    }
    edge_name <- gsub('"', "", parts[1], fixed = TRUE)
    if (edge_name %in% edge_ids) {
      hit_lines <- c(hit_lines, line)
    }
    if (length(hit_lines) == length(edge_ids)) {
      break
    }
  }

  fread(text = c(header_line, hit_lines))
}

cat("Step 1: load score tables and case metadata\n")
score_wide <- read.csv(
  file.path(out_dir, "module_scores_wide.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
edge_membership <- read.csv(
  file.path(out_dir, "module_edge_membership.csv"),
  stringsAsFactors = FALSE
)

case_row <- score_wide %>%
  filter(Sample == case_sample)

if (nrow(case_row) != 1) {
  stop("Case sample was not found in module_scores_wide.csv")
}

cat("Step 2: build cohort scatter plot with highlighted sample\n")
scatter_df <- score_wide %>%
  filter(Response %in% c("CR", "PR", "SD", "PD")) %>%
  mutate(
    is_case = ifelse(Sample == case_sample, "Case sample", "Other samples")
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = expr_Checkpoint_neighborhood,
    y = int_Checkpoint_neighborhood,
    color = Response
  )
) +
  geom_vline(
    xintercept = median(scatter_df$expr_Checkpoint_neighborhood, na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey65"
  ) +
  geom_hline(
    yintercept = median(scatter_df$int_Checkpoint_neighborhood, na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey65"
  ) +
  geom_point(alpha = 0.65, size = 2) +
  geom_point(
    data = case_row,
    aes(
      x = expr_Checkpoint_neighborhood,
      y = int_Checkpoint_neighborhood
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 4.5,
    stroke = 1.2,
    fill = "gold",
    color = "black"
  ) +
  geom_text(
    data = case_row,
    aes(
      x = expr_Checkpoint_neighborhood,
      y = int_Checkpoint_neighborhood,
      label = "SAM2eb07dedf07f"
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = -0.6,
    size = 3.1
  ) +
  scale_color_manual(
    values = c(CR = "#2C7BB6", PR = "#ABD9E9", SD = "#FDAE61", PD = "#D7191C")
  ) +
  labs(
    title = "Case sample in the checkpoint expression-interaction landscape",
    x = "Checkpoint expression score",
    y = "Checkpoint interaction score",
    color = "Response"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

cat("Step 3: extract checkpoint-neighborhood edge weights for the case sample\n")
checkpoint_edges <- edge_membership %>%
  filter(module == "Checkpoint_neighborhood")

matrix_path <- file.path(base_dir, "immune/1Matrix/MOSSN_noPrior_Matrix.csv")
all_link_ids <- fread(matrix_path, select = "Interaction")[[1]]
all_parts <- strsplit(all_link_ids, "_", fixed = TRUE)
all_ok <- lengths(all_parts) >= 2
all_edge_df <- data.frame(
  edge = all_link_ids[all_ok],
  gene1 = vapply(all_parts[all_ok], function(z) z[1], character(1)),
  gene2 = vapply(all_parts[all_ok], function(z) z[2], character(1)),
  stringsAsFactors = FALSE
)

checkpoint_core_extra_edges <- all_edge_df %>%
  filter(
    gene1 %in% checkpoint_core_genes,
    gene2 %in% checkpoint_core_genes
  ) %>%
  mutate(
    module = "Checkpoint_core_augmented",
    edge_type = "checkpoint_core"
  )

case_edge_membership <- checkpoint_core_extra_edges

edge_dt <- extract_edge_matrix_subset(
  matrix_path = matrix_path,
  edge_ids = case_edge_membership$edge
)

edge_weight_df <- data.frame(
  edge = edge_dt$Interaction,
  weight = as.numeric(edge_dt[[case_sample]]),
  stringsAsFactors = FALSE
) %>%
  inner_join(case_edge_membership, by = "edge") %>%
  filter(weight > 0) %>%
  arrange(desc(weight))

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

edge_weight_df <- edge_weight_df %>%
  filter(edge %in% keep_edges)

cat("Step 4: extract case-sample gene expression for network nodes\n")
expr_dt <- fread(file.path(base_dir, "immune/IM210/data/exprSet_filtered.csv"))
expr_genes <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- expr_genes

sample_ids <- colnames(expr_mat)
expr_z <- t(apply(expr_mat, 1, function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}))
rownames(expr_z) <- rownames(expr_mat)
colnames(expr_z) <- sample_ids

node_genes <- sort(unique(c(edge_weight_df$gene1, edge_weight_df$gene2)))
expr_map <- data.frame(
  gene_expr = rownames(expr_z),
  gene_norm = normalize_gene_symbol(rownames(expr_z)),
  stringsAsFactors = FALSE
) %>%
  distinct(gene_norm, .keep_all = TRUE)

node_df <- data.frame(
  gene = node_genes,
  gene_norm = normalize_gene_symbol(node_genes),
  stringsAsFactors = FALSE
) %>%
  left_join(expr_map, by = "gene_norm") %>%
  mutate(
    expr_z = expr_z[gene_expr, case_sample],
    node_role = ifelse(
      gene %in% c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD80", "CD86", "LAG3",
                  "TIGIT", "HAVCR2", "CD28", "ICOS", "ICOSLG", "TNFRSF9",
                  "TNFRSF4", "CD40", "CD40LG", "BTLA", "VSIR", "ENTPD1",
                  "LAIR1", "CD27", "CD70"),
      "Checkpoint gene",
      "Immune partner"
    )
  )

graph_obj <- graph_from_data_frame(
  d = edge_weight_df[, c("gene1", "gene2", "weight")],
  vertices = node_df,
  directed = FALSE
)

set.seed(123)
p_network <- ggraph(graph_obj, layout = "fr") +
  geom_edge_link(
    aes(width = weight, color = weight),
    alpha = 0.8,
    show.legend = TRUE
  ) +
  scale_edge_width_continuous(range = c(0.5, 3)) +
  scale_edge_color_gradient(low = "#FEE0D2", high = "#CB181D") +
  geom_node_point(
    aes(fill = expr_z, shape = node_role),
    size = 7,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  scale_shape_manual(values = c(`Checkpoint gene` = 21, `Immune partner` = 22)) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3.5,
    fontface = "bold"
  ) +
  labs(
    title = "Top checkpoint-neighborhood edges in the case sample",
    edge_color = "MOSSN weight",
    edge_width = "MOSSN weight",
    fill = "Gene expression z-score",
    shape = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

node_attr <- node_df[, c("gene", "expr_z", "node_role")]
cytoscape_df <- edge_weight_df[, c("gene1", "gene2", "weight", "edge_type")] %>%
  left_join(node_attr, by = c("gene1" = "gene")) %>%
  rename(gene1_expr_z = expr_z, gene1_role = node_role) %>%
  left_join(node_attr, by = c("gene2" = "gene")) %>%
  rename(gene2_expr_z = expr_z, gene2_role = node_role)

cytoscape_df <- cytoscape_df %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

nodes1 <- data.frame(gene = node_df$gene, expr_z = node_df$expr_z,
                     role = node_df$node_role)
nodes2 <- data.frame(gene = node_df$gene, expr_z = node_df$expr_z,
                     role = node_df$node_role)
node_out <- unique(rbind(nodes1, nodes2))

write.csv(cytoscape_df, file.path(out_dir, "cytoscape_network.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(node_out, file.path(out_dir, "cytoscape_node_attr.csv"),
          row.names = FALSE, quote = FALSE)
cat("Finished case-sample visualization.\n")
