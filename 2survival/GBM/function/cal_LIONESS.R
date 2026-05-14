rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- "GBM"

available_datasets <- c("CGGA301", "CGGA325", "CGGA693",
                          "GSE4412", "GSE13041", "GSE16011",
                          "GSE72951", "GSE74187", "GSE83300",
                          "TCGAGBM")


source('/proj/c.zihao/work1/function/LIONESS.R')

for (disease_name in available_datasets) {
    print(disease_name)

    save_path <- paste0("/proj/c.zihao/work1/2survival/", researchAim, '/', disease_name, "/LIONESS/")
    ppiFile   <- "/proj/c.zihao/work1/1NT/1data/string/links.csv"
    exprSetFile <- paste0("/proj/c.zihao/work1/2survival/", researchAim,  '/', disease_name, "/data/", "exprSet_filtered.csv")

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

    resultDF <- LIONcal(exprSetFile, ppiFile)

    setwd(save_path)
    write.csv(resultDF, file = "result.csv", row.names = FALSE)
}
