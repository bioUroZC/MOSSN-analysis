rm(list=ls())

cancertypes <- c("BLCA", "LIHC", "LUAD", "SARC", "STAD")

source_dir <- "/proj/c.zihao/work1/04multi/1prepare/05EXPfiles"
output_dir <- "/proj/c.zihao/work1/04multi/1prepare/05EXPout"
link_file <- "/proj/c.zihao/work1/01data/06string/links.csv"
links <- read.csv(link_file, row.names = 1)
ppi_genes <- unique(c(links$protein1, links$protein2))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (tumor in cancertypes) {
  source_file <- file.path(source_dir, tumor, "exprSet.csv")

  if (!file.exists(source_file)) {
    stop(sprintf("Expression file not found: %s", source_file))
  }

  target_file <- file.path(
    output_dir,
    sprintf("EXP_%s.csv", tumor)
  )

  exp_data <- read.csv(source_file, row.names = 1, check.names = FALSE)
  exp_data <- exp_data[rownames(exp_data) %in% ppi_genes, , drop = FALSE]
  exp_data <- log2(exp_data + 1)
  exp_data <- round(exp_data, 5)
  write.csv(exp_data, target_file)
  message(sprintf("Saved rounded expression data -> %s", target_file))
}
