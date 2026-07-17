rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
set.seed(42)

researchAim <- '4net/2intact'

available_datasets <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                        "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

source(paste0(PROJ_ROOT, '/function/LIONESS.R'))

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- paste0(paste0(PROJ_ROOT, "/1NT/"), researchAim, "/LIONESS/", disease_name)
    ppiFile   <- paste0(PROJ_ROOT, "/1NT/4net/2intact/intact_link.csv")
    exprSetFile <- paste0(paste0(PROJ_ROOT, "/1NT/1data/exprset/"), disease_name, "_exprSet_filtered.csv")

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
