rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
set.seed(42)

researchAim <- "6robust/50"
dataset_name <- "LUAD"

source(paste0(PROJ_ROOT, '/function/LIONESS.R'))

print(dataset_name)

save_path <- paste0(paste0(PROJ_ROOT, "/1NT/"), researchAim, "/LIONESS/", dataset_name)
ppiFile   <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")
exprSetFile <- paste0(paste0(PROJ_ROOT, "/1NT/"), researchAim, "/data/LUAD_exprSet_half.csv")

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)

resultDF <- LIONcal(exprSetFile, ppiFile)

setwd(save_path)
write.csv(resultDF, file = "result.csv", row.names = FALSE)
