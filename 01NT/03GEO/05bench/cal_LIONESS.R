rm(list = ls())
library(dplyr)
set.seed(42)

base_dir <- "/proj/c.zihao/work1/01NT/03GEO/Output"

available_datasets <- c("LUAD", "eLUAD", "KIRC", "PRAD")
dataset_dirs <- c(LUAD = "01LUAD", eLUAD = "02eLUAD", KIRC = "03KIRC", PRAD = "04PRAD")

source('/proj/c.zihao/work1/function/benchmark/LIONESS.R')

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- file.path(base_dir, "LIONESS", disease_name)
    ppiFile   <- "/proj/c.zihao/work1/01data/06string/links.csv"
    exprSetFile <- file.path(
        "/proj/c.zihao/work1/01NT/03GEO",
        dataset_dirs[[disease_name]],
        paste0(disease_name, "_exprSet_filtered.csv")
    )

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
