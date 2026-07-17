rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(igraph)
})

set.seed(1)

base_dir <- paste0(PROJ_ROOT, "/1NT/9analysis")
setwd(base_dir)

atlas_file <- file.path(base_dir, "atlas_all.csv")
recurrent_file <- file.path(base_dir, "universal_recurrent_links.csv")
out_dir <- file.path(base_dir, "module_results")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cancers <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
             "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

atlas_all <- read.csv(
  atlas_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

split_link <- function(link_vec) {
  parts <- stringr::str_split_fixed(link_vec, "_", 2)
  data.frame(
    gene1 = parts[, 1],
    gene2 = parts[, 2],
    stringsAsFactors = FALSE
  )
}

edge_weight_rescale <- function(x) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  if (all(is.na(x))) {
    return(rep(1, length(x)))
  }
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) < 1e-8) {
    return(rep(1, length(x)))
  }
  1 + 4 * (x - rng[1]) / diff(rng)
}

run_module_detection <- function(link_file, direction_label) {
  class_label <- if (direction_label == "gain") {
    "recurrently_gained"
  } else {
    "recurrently_lost"
  }
  
  links_df <- read.csv(link_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    dplyr::filter(recurrent_class == class_label)
  if (nrow(links_df) == 0) {
    return(NULL)
  }
  
  genes_df <- split_link(links_df$link)
  edges <- dplyr::bind_cols(links_df, genes_df) %>%
    dplyr::filter(
      gene1 != "",
      gene2 != "",
      !is.na(gene1),
      !is.na(gene2),
      gene1 != gene2
    ) %>%
    dplyr::mutate(
      weight = recurrent_count * abs(median_delta_across_cancers)
    ) %>%
    dplyr::arrange(dplyr::desc(weight), dplyr::desc(recurrent_count))
  
  graph_obj <- igraph::graph_from_data_frame(
    d = edges[, c("gene1", "gene2", "weight", "recurrent_count", "median_delta_across_cancers", "link"), drop = FALSE],
    directed = FALSE
  )
  
  graph_obj <- igraph::simplify(
    graph_obj,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(
      weight = "max",
      recurrent_count = "max",
      median_delta_across_cancers = "mean",
      link = "first"
    )
  )
  
  if (igraph::ecount(graph_obj) == 0) {
    return(NULL)
  }
  
  comm <- igraph::cluster_louvain(graph_obj, weights = igraph::E(graph_obj)$weight)
  membership_vec <- igraph::membership(comm)
  membership_tbl <- data.frame(
    gene = names(membership_vec),
    community_id = as.integer(membership_vec),
    stringsAsFactors = FALSE
  )
  
  membership_tbl <- membership_tbl %>%
    dplyr::count(community_id, name = "n_genes") %>%
    dplyr::arrange(dplyr::desc(n_genes), community_id) %>%
    dplyr::mutate(
      module_rank = dplyr::row_number(),
      module_id = paste0(direction_label, "_M", sprintf("%02d", module_rank))
    ) %>%
    dplyr::select(community_id, module_id, module_rank, n_genes) %>%
    dplyr::right_join(
      data.frame(
        gene = names(membership_vec),
        community_id = as.integer(membership_vec),
        stringsAsFactors = FALSE
      ),
      by = "community_id"
    )
  
  gene_tbl <- membership_tbl %>%
    dplyr::group_by(module_id, module_rank, gene) %>%
    dplyr::summarise(.groups = "drop")
  
  edge_tbl <- edges %>%
    dplyr::left_join(
      membership_tbl %>% dplyr::select(gene, module_id, module_rank),
      by = c("gene1" = "gene")
    ) %>%
    dplyr::rename(module_id_1 = module_id, module_rank_1 = module_rank) %>%
    dplyr::left_join(
      membership_tbl %>% dplyr::select(gene, module_id, module_rank),
      by = c("gene2" = "gene")
    ) %>%
    dplyr::rename(module_id_2 = module_id, module_rank_2 = module_rank) %>%
    dplyr::filter(module_id_1 == module_id_2) %>%
    dplyr::mutate(
      module_id = module_id_1,
      module_rank = module_rank_1,
      direction = direction_label
    ) %>%
    dplyr::select(
      module_id, module_rank, direction, link, gene1, gene2,
      recurrent_count, median_delta_across_cancers, weight,
      n_gain, n_loss, consistency, cancers_gain, cancers_loss
    )
  
  module_summary <- edge_tbl %>%
    dplyr::group_by(module_id, module_rank, direction) %>%
    dplyr::summarise(
      n_edges = dplyr::n(),
      n_genes = dplyr::n_distinct(c(gene1, gene2)),
      median_recurrence = stats::median(recurrent_count, na.rm = TRUE),
      median_delta = stats::median(median_delta_across_cancers, na.rm = TRUE),
      total_weight = sum(weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  gene_stats <- gene_tbl %>%
    dplyr::left_join(
      edge_tbl %>%
        dplyr::select(module_id, gene1, gene2, weight) %>%
        tidyr::pivot_longer(
          cols = c(gene1, gene2),
          names_to = "which_end",
          values_to = "gene"
        ) %>%
        dplyr::group_by(module_id, gene) %>%
        dplyr::summarise(
          module_degree = dplyr::n(),
          module_strength = sum(weight, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("module_id", "gene")
    )
  
  top_hubs <- gene_stats %>%
    dplyr::arrange(module_id, dplyr::desc(module_degree), dplyr::desc(module_strength), gene) %>%
    dplyr::group_by(module_id) %>%
    dplyr::summarise(
      top_hubs = paste(head(gene, 5), collapse = ", "),
      .groups = "drop"
    )
  
  module_summary <- module_summary %>%
    dplyr::left_join(top_hubs, by = "module_id") %>%
    dplyr::arrange(module_rank)
  
  message(direction_label, ": ", nrow(module_summary), " modules detected")
  message(direction_label, ": top module = ", module_summary$module_id[1],
          " (", module_summary$n_genes[1], " genes, ",
          module_summary$n_edges[1], " edges)")
  
  list(
    graph = graph_obj,
    edge_tbl = edge_tbl,
    gene_tbl = gene_stats,
    module_summary = module_summary
  )
}

gain_res <- run_module_detection(recurrent_file, "gain")
loss_res <- run_module_detection(recurrent_file, "loss")

combined_edges <- dplyr::bind_rows(
  if (!is.null(gain_res)) gain_res$edge_tbl,
  if (!is.null(loss_res)) loss_res$edge_tbl
)

combined_genes <- dplyr::bind_rows(
  if (!is.null(gain_res)) gain_res$gene_tbl,
  if (!is.null(loss_res)) loss_res$gene_tbl
)

combined_summary <- dplyr::bind_rows(
  if (!is.null(gain_res)) gain_res$module_summary,
  if (!is.null(loss_res)) loss_res$module_summary
)

if (nrow(combined_edges) > 0) {
  utils::write.csv(
    combined_edges,
    file.path(out_dir, "module_edges.csv"),
    row.names = FALSE
  )
}

if (nrow(combined_genes) > 0) {
  utils::write.csv(
    combined_genes,
    file.path(out_dir, "module_genes.csv"),
    row.names = FALSE
  )
}

if (nrow(combined_summary) > 0) {
  utils::write.csv(
    combined_summary,
    file.path(out_dir, "module_summary.csv"),
    row.names = FALSE
  )
}

message("Module detection finished. Results saved in: ", out_dir)
