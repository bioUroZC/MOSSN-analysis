rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- "CRC"

available_datasets <- c("ACICAM", "GSE12945", "GSE17536",
                          "GSE17537", "GSE28722", "GSE29621",
                          "GSE39582", "GSE41258", "TCGACRC")


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
