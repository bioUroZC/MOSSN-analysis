rm(list = ls())
library(dplyr)
set.seed(42)

base_dir <- "/proj/c.zihao/work1/01NT/01string/02benchmark"

available_datasets <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                        "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

source('/proj/c.zihao/work1/function/benchmark/LIONESS.R')

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- file.path(base_dir, "LIONESS", disease_name)
    ppiFile   <- "/proj/c.zihao/work1/01data/06string/links.csv"
    exprSetFile <- file.path(
        "/proj/c.zihao/work1/01data/exprset",
        paste0(disease_name, "_exprSet_filtered.csv")
    )

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
