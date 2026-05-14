rm(list = ls())

library(data.table)
library(dplyr)

base_dir <- "/proj/c.zihao/work1/3drugs"
out_dir <- file.path(base_dir, "immune/IM210/case/results")

module_order <- c(
  "IFNG_response",
  "Antigen_presentation",
  "Cytotoxicity",
  "Checkpoint_neighborhood"
)

row_zscore <- function(mat) {
  t(apply(mat, 1, function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
      return(rep(0, length(x)))
    }
    as.numeric(scale(x))
  }))
}

cat("Step 1: load clinical cohort\n")
follow_up <- read.csv(
  file.path(base_dir, "immune/IM210/data/IMvigor210_FollowUp.csv"),
  header = TRUE,
  row.names = 1
)
follow_up$Sample <- rownames(follow_up)
follow_up <- subset(
  follow_up,
  select = c(
    "Sample",
    "Best.Confirmed.Overall.Response",
    "binaryResponse",
    "Immune.phenotype",
    "Tissue",
    "os",
    "censOS"
  )
)
names(follow_up) <- c("Sample", "Response", "Binary", "Phenotype", "Tissue", "Time", "OS")

cohort_df <- follow_up %>%
  filter(
    Tissue == "bladder",
    !is.na(Response),
    Response %in% c("CR", "PR", "SD", "PD")
  ) %>%
  mutate(
    responder = ifelse(Response %in% c("CR", "PR"), "Responder", "Non_responder"),
    response_ordinal = c(CR = 4, PR = 3, SD = 2, PD = 1)[Response]
  )

sample_ids <- cohort_df$Sample
cat("Bladder samples retained:", length(sample_ids), "\n")

cat("Step 2: load predefined modules\n")
gene_membership <- read.csv(
  file.path(out_dir, "module_gene_membership.csv"),
  stringsAsFactors = FALSE
)
edge_membership <- read.csv(
  file.path(out_dir, "module_edge_membership.csv"),
  stringsAsFactors = FALSE
)

cat("Step 3: compute expression-level module scores\n")
expr_dt <- fread(file.path(base_dir, "immune/IM210/data/exprSet_filtered.csv"))
expr_genes <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- expr_genes

expr_mat <- expr_mat[, colnames(expr_mat) %in% sample_ids, drop = FALSE]
expr_mat <- expr_mat[, sample_ids, drop = FALSE]
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- expr_genes
colnames(expr_mat) <- sample_ids

expr_z <- row_zscore(expr_mat)
rownames(expr_z) <- rownames(expr_mat)
colnames(expr_z) <- colnames(expr_mat)

expr_gene_sets <- gene_membership %>%
  filter(module %in% module_order, in_expression) %>%
  group_by(module) %>%
  summarise(genes = list(unique(gene_expr)), .groups = "drop")

expr_score_list <- lapply(seq_len(nrow(expr_gene_sets)), function(i) {
  module_name <- expr_gene_sets$module[i]
  genes <- expr_gene_sets$genes[[i]]
  score <- colMeans(expr_z[genes, , drop = FALSE], na.rm = TRUE)
  data.frame(
    Sample = names(score),
    module = module_name,
    expression_score = as.numeric(score),
    stringsAsFactors = FALSE
  )
})

expr_scores_long <- bind_rows(expr_score_list)

cat("Step 4: extract module edges from MOSSN matrix\n")
edge_keep <- unique(edge_membership$edge)
edge_file <- tempfile(fileext = ".txt")
writeLines(edge_keep, edge_file)

matrix_path <- file.path(base_dir, "immune/1Matrix/MOSSN_noPrior_Matrix.csv")
awk_cmd <- sprintf(
  "awk -F',' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} {key=$1; gsub(/\"/, \"\", key); if (key in keep) print}' %s %s",
  shQuote(edge_file),
  shQuote(matrix_path)
)

link_dt <- fread(cmd = awk_cmd)
unlink(edge_file)

link_ids <- link_dt$Interaction
link_mat <- as.matrix(link_dt[, -1, with = FALSE])
rownames(link_mat) <- link_ids

link_mat <- link_mat[, colnames(link_mat) %in% sample_ids, drop = FALSE]
link_mat <- link_mat[, sample_ids, drop = FALSE]
link_mat <- apply(link_mat, 2, as.numeric)
rownames(link_mat) <- link_ids
colnames(link_mat) <- sample_ids

keep_rows <- rowSums(link_mat != 0, na.rm = TRUE) > 0
link_mat <- link_mat[keep_rows, , drop = FALSE]
cat("Module edges extracted (non-zero in IM210):", nrow(link_mat), "\n")

cat("Step 5: compute interaction-level module scores\n")
edge_sets <- edge_membership %>%
  group_by(module) %>%
  summarise(edges = list(unique(edge)), .groups = "drop")

interaction_score_list <- lapply(seq_len(nrow(edge_sets)), function(i) {
  module_name <- edge_sets$module[i]
  edges <- edge_sets$edges[[i]]
  score <- colMeans(link_mat[edges, , drop = FALSE], na.rm = TRUE)
  data.frame(
    Sample = names(score),
    module = module_name,
    interaction_score = as.numeric(score),
    stringsAsFactors = FALSE
  )
})

interaction_scores_long <- bind_rows(interaction_score_list)

cat("Step 6: combine scores and build wide table\n")
score_long <- expr_scores_long %>%
  inner_join(interaction_scores_long, by = c("Sample", "module"))

expr_wide <- reshape(
  score_long[, c("Sample", "module", "expression_score")],
  idvar = "Sample",
  timevar = "module",
  direction = "wide"
)
int_wide <- reshape(
  score_long[, c("Sample", "module", "interaction_score")],
  idvar = "Sample",
  timevar = "module",
  direction = "wide"
)

colnames(expr_wide) <- sub("^expression_score\\.", "expr_", colnames(expr_wide))
colnames(int_wide) <- sub("^interaction_score\\.", "int_", colnames(int_wide))

score_wide <- cohort_df %>%
  left_join(expr_wide, by = "Sample") %>%
  left_join(int_wide, by = "Sample")

expr_cols <- paste0("expr_", module_order)
int_cols <- paste0("int_", module_order)

score_wide[, expr_cols] <- scale(score_wide[, expr_cols])
score_wide[, int_cols] <- scale(score_wide[, int_cols])

score_wide$combined_expression_score <- rowMeans(score_wide[, expr_cols], na.rm = TRUE)
score_wide$combined_interaction_score <- rowMeans(score_wide[, int_cols], na.rm = TRUE)

write.csv(
  score_long,
  file.path(out_dir, "module_scores_long.csv"),
  row.names = FALSE
)
write.csv(
  score_wide,
  file.path(out_dir, "module_scores_wide.csv"),
  row.names = FALSE
)

cat("Finished module score computation.\n")
cat("Samples scored:", nrow(score_wide), "\n")
print(table(score_wide$Response))
