rm(list = ls())

library(igraph)

base_dir <- "/proj/c.zihao/work1/05atlas"
out_dir <- file.path(base_dir, "module_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)  

recurrent_tbl <- read.csv(file.path(base_dir, "universal_recurrent_links.csv"))

all_edges <- data.frame()
all_genes <- data.frame()
all_summary <- data.frame()

for (direction in c("gained", "lost")) {

  if (direction == "gained") {
    edges <- recurrent_tbl[recurrent_tbl$recurrent_class == "recurrently_gained", ]
  } else {
    edges <- recurrent_tbl[recurrent_tbl$recurrent_class == "recurrently_lost", ]
  }
  edges$gene1 <- sub("_.*$", "", edges$link)
  edges$gene2 <- sub("^[^_]*_", "", edges$link)
  edges <- edges[edges$gene1 != edges$gene2, ]

  edges$weight <- edges$recurrent_count * abs(edges$median_delta_across_cancers)
  edges <- edges[order(-edges$weight, -edges$recurrent_count), ]

  graph_obj <- graph_from_data_frame(
    edges[, c("gene1", "gene2", "weight", "recurrent_count",
              "median_delta_across_cancers", "link")],
    directed = FALSE
  )
  graph_obj <- simplify(
    graph_obj,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "max", recurrent_count = "max",
                          median_delta_across_cancers = "mean", link = "first")
  )
  comm <- cluster_louvain(graph_obj, weights = E(graph_obj)$weight)

  genes <- data.frame(gene = names(membership(comm)),
                      community = as.integer(membership(comm)))

  size_tbl <- as.data.frame(table(community = genes$community))
  size_tbl$community <- as.integer(as.character(size_tbl$community))
  size_tbl <- size_tbl[order(-size_tbl$Freq, size_tbl$community), ]
  size_tbl$module_rank <- seq_len(nrow(size_tbl))
  size_tbl$module_id <- paste0(direction, "_M", sprintf("%02d", size_tbl$module_rank))

  idx <- match(genes$community, size_tbl$community)
  genes$module_id <- size_tbl$module_id[idx]
  genes$module_rank <- size_tbl$module_rank[idx]

  module1 <- genes$module_id[match(edges$gene1, genes$gene)]
  module2 <- genes$module_id[match(edges$gene2, genes$gene)]
  edges <- edges[module1 == module2, ]
  edges$module_id <- genes$module_id[match(edges$gene1, genes$gene)]
  edges$module_rank <- genes$module_rank[match(edges$gene1, genes$gene)]
  edges$direction <- direction
  edges <- edges[, c("module_id", "module_rank", "direction", "link", "gene1", "gene2",
                     "recurrent_count", "median_delta_across_cancers", "weight",
                     "n_gained", "n_lost", "consistency", "cancers_gained", "cancers_lost")]

  genes <- genes[order(genes$module_id, genes$gene, method = "radix"), ]
  genes$module_degree <- NA
  genes$module_strength <- NA
  for (i in seq_len(nrow(genes))) {
    in_module <- edges$module_id == genes$module_id[i]
    touches <- in_module & (edges$gene1 == genes$gene[i] | edges$gene2 == genes$gene[i])
    if (any(touches)) {
      genes$module_degree[i] <- sum(touches)
      genes$module_strength[i] <- sum(edges$weight[touches])
    }
  }
  genes <- genes[, c("module_id", "module_rank", "gene", "module_degree", "module_strength")]

  summary_tbl <- data.frame()
  for (module_id in sort(unique(edges$module_id), method = "radix")) {
    e <- edges[edges$module_id == module_id, ]
    g <- genes[genes$module_id == module_id, ]

    g <- g[order(-g$module_degree, -g$module_strength, g$gene, method = "radix"), ]

    summary_tbl <- rbind(summary_tbl, data.frame(
      module_id = module_id,
      module_rank = e$module_rank[1],
      direction = direction,
      n_edges = nrow(e),
      n_genes = length(unique(c(e$gene1, e$gene2))),
      median_recurrence = median(e$recurrent_count),
      median_delta = median(e$median_delta_across_cancers),
      total_weight = sum(e$weight),
      top_hubs = paste(head(g$gene, 5), collapse = ", ")
    ))
  }
  summary_tbl <- summary_tbl[order(summary_tbl$module_rank), ]

  message(direction, ": ", nrow(summary_tbl), " modules detected")
  message(direction, ": top module = ", summary_tbl$module_id[1],
          " (", summary_tbl$n_genes[1], " genes, ", summary_tbl$n_edges[1], " edges)")

  all_edges <- rbind(all_edges, edges)
  all_genes <- rbind(all_genes, genes)
  all_summary <- rbind(all_summary, summary_tbl)
}

write.csv(all_edges, file.path(out_dir, "module_edges.csv"), row.names = FALSE)
write.csv(all_genes, file.path(out_dir, "module_genes.csv"), row.names = FALSE)
write.csv(all_summary, file.path(out_dir, "module_summary.csv"), row.names = FALSE)
