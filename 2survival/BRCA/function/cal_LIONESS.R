rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- "BRCA"

available_datasets <- c("GSE11121", "GSE12093", "GSE162228",
                          "GSE17705", "GSE20685", "GSE20711",
                          "GSE21653", "GSE22219", "GSE25055",
                          "GSE25065", "GSE42568", "GSE45255",
                          "GSE48390", "GSE61304", "GSE7390",
                          "TCGABRCA")


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
