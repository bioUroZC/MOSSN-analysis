rm(list = ls())

base_dir <- "/proj/c.zihao/work1/06case"
out_dir <- file.path(base_dir, "04results")

clinical_file <- file.path(base_dir, "01data/IMvigor210_FollowUp.csv")
expr_file <- file.path(base_dir, "01data/exprSet_filtered.csv")
matrix_file <- file.path(base_dir, "04results/MOSSN/merged_matrix.csv")

module_order <- c("IFNG_response", "Antigen_presentation",
                  "Cytotoxicity", "Checkpoint_neighborhood")

clinical <- read.csv(clinical_file, row.names = 1)
cohort <- data.frame(
  Sample = rownames(clinical),
  Response = clinical$Best.Confirmed.Overall.Response,
  Binary = clinical$binaryResponse,
  Phenotype = clinical$Immune.phenotype,
  Tissue = clinical$Tissue,
  Time = clinical$os,
  OS = clinical$censOS
)
keep <- cohort$Tissue %in% "bladder" & cohort$Response %in% c("CR", "PR", "SD", "PD")
cohort <- cohort[keep, ]

cohort$responder <- ifelse(cohort$Response %in% c("CR", "PR"), "Responder", "Non_responder")
cohort$response_ordinal <- c(CR = 4, PR = 3, SD = 2, PD = 1)[cohort$Response]

sample_ids <- cohort$Sample
cat("Bladder samples retained:", length(sample_ids), "\n")

gene_membership <- read.csv(file.path(out_dir, "module_gene_membership.csv"))
edge_membership <- read.csv(file.path(out_dir, "module_edge_membership.csv"))

modules <- sort(unique(edge_membership$module))

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr[, sample_ids])

gene_sd <- apply(expr_mat, 1, sd)
expr_z <- (expr_mat - rowMeans(expr_mat)) / gene_sd
expr_z[gene_sd == 0, ] <- 0

link <- data.table::fread(matrix_file, data.table = FALSE)
link <- link[link$Interaction %in% edge_membership$edge, ]
link_mat <- as.matrix(link[, sample_ids])
rownames(link_mat) <- link$Interaction

link_mat <- link_mat[rowSums(link_mat != 0) > 0, ]
cat("Module edges extracted (non-zero in IM210):", nrow(link_mat), "\n")

score_long <- data.frame()
for (module_name in modules) {
  in_module <- gene_membership$module == module_name & gene_membership$in_expression
  genes <- unique(gene_membership$gene_expr[in_module])
  edges <- unique(edge_membership$edge[edge_membership$module == module_name])

  score_long <- rbind(score_long, data.frame(
    Sample = sample_ids,
    module = module_name,
    expression_score = colMeans(expr_z[genes, , drop = FALSE]),
    interaction_score = colMeans(link_mat[edges, , drop = FALSE])
  ))
}
rownames(score_long) <- NULL

score_wide <- cohort
for (module_name in modules) {
  one <- score_long[score_long$module == module_name, ]
  score_wide[[paste0("expr_", module_name)]] <- one$expression_score[match(score_wide$Sample, one$Sample)]
}
for (module_name in modules) {
  one <- score_long[score_long$module == module_name, ]
  score_wide[[paste0("int_", module_name)]] <- one$interaction_score[match(score_wide$Sample, one$Sample)]
}

expr_cols <- paste0("expr_", module_order)
int_cols <- paste0("int_", module_order)
score_wide[, expr_cols] <- scale(score_wide[, expr_cols])
score_wide[, int_cols] <- scale(score_wide[, int_cols])

score_wide$combined_expression_score <- rowMeans(score_wide[, expr_cols])
score_wide$combined_interaction_score <- rowMeans(score_wide[, int_cols])

num_cols <- sapply(score_long, is.numeric)
score_long[num_cols] <- round(score_long[num_cols], 5)
num_cols <- sapply(score_wide, is.numeric)
score_wide[num_cols] <- round(score_wide[num_cols], 5)

write.csv(score_long, file.path(out_dir, "module_scores_long.csv"), row.names = FALSE)
write.csv(score_wide, file.path(out_dir, "module_scores_wide.csv"), row.names = FALSE)

cat("Samples scored:", nrow(score_wide), "\n")
print(table(score_wide$Response))
