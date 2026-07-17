# Install the R packages required by the MOSSN analysis workflows.
#
# The analyses reported in the manuscript were run on R 4.3.1.
# Usage:
#     Rscript install_R_packages.R
#
# Packages already present are skipped. Annotation packages are large; the
# Bioconductor step may take a while on a fresh installation.

cran_packages <- c(
  "RColorBrewer", "RCurl", "UpSetR", "caret", "cowplot", "data.table",
  "dplyr", "e1071", "flextable", "gdata", "ggplot2", "ggraph", "glmnet",
  "igraph", "msigdbr", "officer", "openxlsx", "pROC", "patchwork",
  "pheatmap", "poweRlaw", "purrr", "randomForest", "readr", "readxl",
  "rvest", "scales", "stringr", "survival", "survminer", "tibble",
  "tidyr", "timeROC"
)

bioc_packages <- c(
  "AnnotationDbi", "ArrayExpress", "Biobase", "GEOquery", "GenomicRanges",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19", "SummarizedExperiment",
  "TCGAbiolinks", "TxDb.Hsapiens.UCSC.hg19.knownGene", "affy", "annotatr",
  "biomaRt", "clusterProfiler", "hgu133plus2.db", "limma", "org.Hs.eg.db",
  "survcomp"
)

missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

to_install <- missing(cran_packages)
if (length(to_install)) {
  message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN packages already installed.")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

to_install <- missing(bioc_packages)
if (length(to_install)) {
  message("Installing Bioconductor packages: ", paste(to_install, collapse = ", "))
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

still_missing <- missing(c(cran_packages, bioc_packages))
if (length(still_missing)) {
  warning("The following packages could not be installed: ",
          paste(still_missing, collapse = ", "))
} else {
  message("All required R packages are available.")
}
