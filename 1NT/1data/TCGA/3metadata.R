# ====================================================================================

# ====================================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)


base_path <- '/proj/c.zihao/work1/1NT/1data/TCGA'

setwd(base_path)

tcga_files <- list.files(
  path = base_path,
  pattern = "TCGA.*\\.csv$",
  full.names = TRUE
)

print(tcga_files)


metadata_list <- list()

for (f in tcga_files) {
  print(f)
  dat <- read.csv(f, header = TRUE, row.names = 1, check.names = FALSE)

  tumor_type <- sub(".*TCGA-([A-Z]+)\\.csv$", "\\1", f)

  meta_tmp <- data.frame(
    Sample = colnames(dat),
    Type    = tumor_type,
    stringsAsFactors = FALSE
  )

  metadata_list[[tumor_type]] <- meta_tmp
}

metadata <- bind_rows(metadata_list)

print(head(metadata))
print(table(metadata$Type))

metadata <- subset(metadata, metadata$Type != "COAD")
metadata <- subset(metadata, metadata$Type != "READ")

print(table(metadata$Type))
setwd(base_path)
write.csv(metadata, file='metadata.csv')
