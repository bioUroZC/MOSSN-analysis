rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- '2intact'

available_datasets <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                        "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

source('/proj/c.zihao/work1/function/LIONESS.R')

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, "/LIONESS/", disease_name)
    ppiFile   <- "/proj/c.zihao/work1/1NT/2intact/intact_link.csv"
    exprSetFile <- paste0("/proj/c.zihao/work1/1NT/1data/exprset/", disease_name, "_exprSet_filtered.csv")

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
