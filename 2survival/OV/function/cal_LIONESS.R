rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- "OV"

available_datasets <- c("GSE102073", "GSE13876", "GSE140082",
                          "GSE17260", "GSE18520", "GSE23554",
                          "GSE26193", "GSE26712", "GSE30161",
                          "GSE31245", "GSE32062", "GSE51088",
                          "GSE53963", "GSE63885", "GSE73614",
                          "GSE8842", "GSE9891", "MTAB386",
                          "TCGAOV")


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
