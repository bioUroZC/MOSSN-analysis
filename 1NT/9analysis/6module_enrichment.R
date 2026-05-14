rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
module_dir <- file.path(base_dir, "module_results")
out_dir <- module_dir

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

module_summary <- read.csv(
  file.path(module_dir, "module_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

module_genes <- read.csv(
  file.path(module_dir, "module_genes.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
module_genes <- module_genes %>%
  dplyr::distinct(module_id, gene, .keep_all = TRUE)

universe_symbols <- sort(unique(module_genes$gene))
universe_map <- clusterProfiler::bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
universe_entrez <- sort(unique(universe_map$ENTREZID))

hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::distinct() %>%
  dplyr::rename(term = gs_name, ENTREZID = ncbi_gene)

clean_term <- function(x) {
  x %>%
    stringr::str_replace("^HALLMARK_", "") %>%
    stringr::str_replace("^KEGG_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_title()
}

pick_top_term <- function(df, source_name) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      source = source_name,
      term = NA_character_,
      fdr = NA_real_,
      count = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  df %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count), Description) %>%
    dplyr::slice(1) %>%
    dplyr::transmute(
      source = source_name,
      term = clean_term(Description),
      fdr = p.adjust,
      count = Count
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

fmt_empty <- function(module_row, reason_text) {
  data.frame(
    module_id = module_row$module_id,
    module_rank = module_row$module_rank,
    direction = module_row$direction,
    n_edges = module_row$n_edges,
    n_genes = module_row$n_genes,
    top_hubs = module_row$top_hubs,
    hallmark_top = NA_character_,
    hallmark_fdr = NA_real_,
    hallmark_count = NA_real_,
    kegg_top = NA_character_,
    kegg_fdr = NA_real_,
    kegg_count = NA_real_,
    annotation_status = reason_text,
    best_source_hk = "Unannotated",
    best_term_hk = "Unannotated",
    concise_label = "Unannotated",
    stringsAsFactors = FALSE
  )
}

all_term_rows <- list()
annotation_rows <- list()

for (i in seq_len(nrow(module_summary))) {
  module_row <- module_summary[i, , drop = FALSE]
  module_id <- module_row$module_id
  
  genes <- module_genes %>%
    dplyr::filter(module_id == !!module_id) %>%
    dplyr::pull(gene) %>%
    unique()
  
  if (length(genes) < 5) {
    annotation_rows[[length(annotation_rows) + 1]] <- fmt_empty(module_row, "insufficient_genes")
    next
  }
  
  gene_map <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  
  if (is.null(gene_map) || nrow(gene_map) < 5) {
    annotation_rows[[length(annotation_rows) + 1]] <- fmt_empty(module_row, "insufficient_mapping")
    next
  }
  
  gene_entrez <- sort(unique(gene_map$ENTREZID))
  
  eh <- clusterProfiler::enricher(
    gene = gene_entrez,
    TERM2GENE = hallmark_df,
    universe = universe_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  eh_df <- as.data.frame(eh)
  
  ekegg <- clusterProfiler::enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    universe = universe_entrez,
    pvalueCutoff = 0.05
  )
  ekegg_df <- as.data.frame(ekegg)
  
  module_terms <- dplyr::bind_rows(
    if (nrow(eh_df) > 0) {
      eh_df %>%
        dplyr::mutate(module_id = module_id, source = "Hallmark") %>%
        dplyr::select(module_id, source, Description, p.adjust, Count)
    },
    if (nrow(ekegg_df) > 0) {
      ekegg_df %>%
        dplyr::mutate(module_id = module_id, source = "KEGG") %>%
        dplyr::select(module_id, source, Description, p.adjust, Count)
    }
  )
  
  if (nrow(module_terms) > 0) {
    all_term_rows[[length(all_term_rows) + 1]] <- module_terms %>%
      dplyr::mutate(clean_description = clean_term(Description))
  }
  
  top_h <- pick_top_term(eh_df, "Hallmark")
  top_k <- pick_top_term(ekegg_df, "KEGG")
  
  best_tbl <- dplyr::bind_rows(top_h, top_k) %>%
    dplyr::filter(!is.na(term)) %>%
    dplyr::arrange(fdr, dplyr::desc(count), term)
  
  if (nrow(best_tbl) == 0) {
    annotation_rows[[length(annotation_rows) + 1]] <- fmt_empty(module_row, "no_significant_terms")
    next
  }
  
  best <- best_tbl[1, , drop = FALSE]
  
  annotation_rows[[length(annotation_rows) + 1]] <- data.frame(
    module_id = module_row$module_id,
    module_rank = module_row$module_rank,
    direction = module_row$direction,
    n_edges = module_row$n_edges,
    n_genes = module_row$n_genes,
    top_hubs = module_row$top_hubs,
    hallmark_top = top_h$term,
    hallmark_fdr = top_h$fdr,
    hallmark_count = top_h$count,
    kegg_top = top_k$term,
    kegg_fdr = top_k$fdr,
    kegg_count = top_k$count,
    annotation_status = "ok",
    best_source_hk = best$source,
    best_term_hk = best$term,
    concise_label = best$term,
    stringsAsFactors = FALSE
  )
}

annotation_tbl <- dplyr::bind_rows(annotation_rows) %>%
  dplyr::arrange(direction, module_rank)

term_tbl <- if (length(all_term_rows) > 0) {
  dplyr::bind_rows(all_term_rows) %>%
    dplyr::arrange(module_id, p.adjust, dplyr::desc(Count), clean_description)
} else {
  data.frame()
}

utils::write.csv(
  annotation_tbl,
  file.path(out_dir, "module_annotation_table.csv"),
  row.names = FALSE
)

utils::write.csv(
  term_tbl,
  file.path(out_dir, "module_enrichment_all_terms.csv"),
  row.names = FALSE
)

plot_df <- annotation_tbl %>%
  dplyr::filter(annotation_status == "ok") %>%
  dplyr::mutate(
    plot_fdr = dplyr::case_when(
      !is.na(hallmark_fdr) & best_source_hk == "Hallmark" ~ hallmark_fdr,
      !is.na(kegg_fdr) & best_source_hk == "KEGG" ~ kegg_fdr,
      TRUE ~ 1
    ),
    score = -log10(pmax(plot_fdr, 1e-300))
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    module_label = paste0(module_id, " | ", best_source_hk, " | ", concise_label),
    module_label = factor(module_label, levels = rev(unique(module_label)))
  )

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = score, y = module_label, fill = direction)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Annotation strength (-log10 FDR)",
      y = "Top module annotations",
      fill = "Direction",
      title = "Top module annotations (Hallmark/KEGG only)"
    )
  
  pdf(file.path(out_dir, "4top_module_annotations.pdf"), width = 10, height = 7)
  print(p)
  dev.off()
}

message("Annotated modules: ", sum(annotation_tbl$annotation_status == "ok"))
message("Saved: ", file.path(out_dir, "module_annotation_table.csv"))
