rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
set.seed(42)

researchAim <- "STAD"

available_datasets <- c("GSE15459", "GSE26253",
                        "GSE26899", "GSE26901", "GSE29272", "GSE57303", "GSE62254",
                       "GSE84437", "TCGASTAD")


source(paste0(PROJ_ROOT, '/function/LIONESS.R'))

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- paste0(paste0(PROJ_ROOT, "/2survival/"), researchAim, '/', disease_name, "/LIONESS/")
    ppiFile   <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")
    exprSetFile <- paste0(paste0(PROJ_ROOT, "/2survival/"), researchAim,  '/', disease_name, "/data/", "exprSet_filtered.csv")

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
