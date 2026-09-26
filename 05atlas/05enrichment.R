rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

base_dir <- "/proj/c.zihao/work1/05atlas"
module_dir <- file.path(base_dir, "module_results")

module_summary <- read.csv(file.path(module_dir, "module_summary.csv"))
module_genes <- read.csv(file.path(module_dir, "module_genes.csv"))

universe_symbols <- sort(unique(module_genes$gene))
universe_map <- bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
universe_entrez <- sort(unique(universe_map$ENTREZID))

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_df <- unique(data.frame(term = hallmark$gs_name, ENTREZID = hallmark$ncbi_gene))

clean_term <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_title(stringr::str_squish(x))
}

annotation_tbl <- data.frame()
term_tbl <- data.frame()

for (i in seq_len(nrow(module_summary))) {
  module_id <- module_summary$module_id[i]
  genes <- unique(module_genes$gene[module_genes$module_id == module_id])

  status <- "insufficient_genes"
  top_h <- data.frame(clean_description = NA, p.adjust = NA, Count = NA)
  top_k <- data.frame(clean_description = NA, p.adjust = NA, Count = NA)
  best_source <- "Unannotated"
  best_term <- "Unannotated"

  if (length(genes) >= 5) {
    gene_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_entrez <- sort(unique(gene_map$ENTREZID))

    hallmark_res <- as.data.frame(enricher(
      gene = gene_entrez, TERM2GENE = hallmark_df, universe = universe_entrez,
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ))
    kegg_res <- as.data.frame(enrichKEGG(
      gene = gene_entrez, organism = "hsa", universe = universe_entrez,
      pvalueCutoff = 0.05
    ))

    cols <- c("Description", "p.adjust", "Count")
    res <- rbind(
      data.frame(module_id = rep(module_id, nrow(hallmark_res)),
                 source = rep("Hallmark", nrow(hallmark_res)),
                 hallmark_res[, cols]),
      data.frame(module_id = rep(module_id, nrow(kegg_res)),
                 source = rep("KEGG", nrow(kegg_res)),
                 kegg_res[, cols])
    )
    res$clean_description <- clean_term(res$Description)
    term_tbl <- rbind(term_tbl, res)

    res <- res[order(res$p.adjust, -res$Count, res$Description, method = "radix"), ]
    top_h <- res[res$source == "Hallmark", ][1, ]
    top_k <- res[res$source == "KEGG", ][1, ]

    best <- rbind(top_h, top_k)
    best <- best[order(best$p.adjust, -best$Count, best$clean_description, method = "radix"), ][1, ]

    status <- ifelse(is.na(best$source), "no_significant_terms", "ok")
    best_source <- ifelse(is.na(best$source), "Unannotated", best$source)
    best_term <- ifelse(is.na(best$source), "Unannotated", best$clean_description)
  }

  annotation_tbl <- rbind(annotation_tbl, data.frame(
    module_id = module_id,
    module_rank = module_summary$module_rank[i],
    direction = module_summary$direction[i],
    n_edges = module_summary$n_edges[i],
    n_genes = module_summary$n_genes[i],
    top_hubs = module_summary$top_hubs[i],
    hallmark_top = top_h$clean_description,
    hallmark_fdr = top_h$p.adjust,
    hallmark_count = top_h$Count,
    kegg_top = top_k$clean_description,
    kegg_fdr = top_k$p.adjust,
    kegg_count = top_k$Count,
    annotation_status = status,
    best_source_hk = best_source,
    best_term_hk = best_term,
    concise_label = best_term
  ))
}

annotation_tbl <- annotation_tbl[order(annotation_tbl$direction, annotation_tbl$module_rank,
                                       method = "radix"), ]

term_tbl <- term_tbl[order(term_tbl$module_id, term_tbl$p.adjust, -term_tbl$Count,
                           term_tbl$clean_description, method = "radix"), ]

write.csv(annotation_tbl, file.path(module_dir, "module_annotation_table.csv"), row.names = FALSE)
write.csv(term_tbl, file.path(module_dir, "module_enrichment_all_terms.csv"), row.names = FALSE)

message("Annotated modules: ", sum(annotation_tbl$annotation_status == "ok"))
