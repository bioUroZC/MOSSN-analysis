
#===================================================

#===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

#===================================================

#===================================================

# URL of the GitHub folder
# github_url <- "https://github.com/camlab-bioml/cancersea/tree/master/data"


setwd(paste0(PROJ_ROOT, "/1NT/1data/Markers"))

rda_files <- list.files(pattern = "\\.rda$", full.names = TRUE)
data_list <- list()

for (file in rda_files) {
  env <- new.env()
  load(file, envir = env)
  obj_name <- ls(env)[1]
  data_list[[obj_name]] <- env[[obj_name]]
}

combined_data <- do.call(rbind, data_list)

data <- combined_data

write.csv(data, file = 'SSmarkers.csv')

