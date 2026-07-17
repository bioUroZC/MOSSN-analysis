#=============================================================

#=============================================================

rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

#=============================================================

#=============================================================

cancertypes <- c("ACC", "BLCA", "BRCA", "CESC", 
                 "CRC", "ESCA", "GBM", "HNSC", "KIRC",
                 "LGG", "LIHC", "LUAD", "LUSC", 
                 "PAAD", "PRAD", "SARC", "STAD")

for (i in 1:length(cancertypes)) {

  setwd(paste0(PROJ_ROOT, '/4coupled/CNVfiles'))
  infile <- paste0('TCGA.', cancertypes[i], '.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz')
  print(cancertypes[i])
  cnv_df <- fread(infile)
  gene_col <- colnames(cnv_df)[1]
  cnv_df <- as.data.frame(cnv_df)

  rownames(cnv_df) <- cnv_df[[gene_col]]
  cnv_df[[gene_col]] <- NULL

  cnv_df[1:5, 1:5]

  print(min(cnv_df))
  print(max(cnv_df))

  cnv_df <- cnv_df + 2

  print(min(cnv_df))
  print(max(cnv_df))

  folder_for_save <- paste0(PROJ_ROOT, "/4coupled/CNVout/")
  dir.create(folder_for_save, showWarnings = FALSE)
  setwd(folder_for_save)

  file_for_save <- paste0("CNV_", cancertypes[i], '.csv')
  write.csv(cnv_df, file = file_for_save)

}
