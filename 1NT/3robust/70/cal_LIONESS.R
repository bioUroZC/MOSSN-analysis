rm(list = ls())
library(dplyr)
set.seed(42)

researchAim <- "3robust/70"
dataset_name <- "LUAD"

source('/proj/c.zihao/work1/function/LIONESS.R')

print(dataset_name)

save_path <- paste0("/proj/c.zihao/work1/1NT/", researchAim, "/LIONESS/", dataset_name)
ppiFile   <- "/proj/c.zihao/work1/1NT/1data/string/links.csv"
exprSetFile <- paste0("/proj/c.zihao/work1/1NT/", researchAim, "/data/LUAD_exprSet_half.csv")

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

resultDF <- LIONcal(exprSetFile, ppiFile)

setwd(save_path)
write.csv(resultDF, file = "result.csv", row.names = FALSE)
