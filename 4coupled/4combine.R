# ======================================================================

# ======================================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(readxl)


# ======================================================================

# ======================================================================

setwd(paste0(PROJ_ROOT, "/2survival"))
survivaldata <- read_excel("mmc1.xlsx")
names(survivaldata)

survivaldata <- subset(survivaldata, select = c("bcr_patient_barcode",
                                                "type","OS" , "OS.time",   "PFI",  "PFI.time" ))

print(table(survivaldata$type))

survivaldata$type[survivaldata$type=="COAD"] <- "CRC"
survivaldata$type[survivaldata$type=="READ"] <- "CRC"


survivaldata1 <- subset(survivaldata, survivaldata$type != "BRCA")
survivaldata2 <- subset(survivaldata, survivaldata$type == "BRCA")


survivaldata1 <- subset(survivaldata1, select = c("bcr_patient_barcode",
                                                "type","OS" , "OS.time"))

survivaldata2 <- subset(survivaldata2, select = c("bcr_patient_barcode",
                                                 "type","PFI",  "PFI.time" ))


names(survivaldata1) <- c("Sample", "Type", "OS", "OSTime")

names(survivaldata2) <- c("Sample", "Type", "OS", "OSTime")

survivaldata <- rbind(survivaldata1, survivaldata2)


cancertypes <- c("ACC", "BLCA", "BRCA", "CESC", 
                 "CRC", "ESCA", "GBM", "HNSC", "KIRC",
                 "LGG", "LIHC", "LUAD", "LUSC", 
                 "PAAD", "PRAD", "SARC", "STAD")


survivaldata <- survivaldata[which(survivaldata$Type %in% cancertypes),]

str(survivaldata)

survivaldata$OS <- as.numeric(survivaldata$OS)

survivaldata$OSTime <- as.numeric(survivaldata$OSTime)

survivaldata$OSTime <- survivaldata$OSTime / 365

survivaldata <- subset(survivaldata, survivaldata$OSTime > 0.1)


survivaldata$Sample <- chartr(old = '.', new = '_', survivaldata$Sample)
survivaldata$Sample <- chartr(old = '-', new = '_', survivaldata$Sample)

print(table(survivaldata$Type))

# ======================================================================

# ======================================================================

save_dir <- paste0(PROJ_ROOT, '/4coupled/files')

exp_dir <- paste0(PROJ_ROOT, '/4coupled/EXPout')

cnv_dir <- paste0(PROJ_ROOT, '/4coupled/CNVout')

met_dir <- paste0(PROJ_ROOT, '/4coupled/METout')


for (cancername in cancertypes) {
  

  setwd(exp_dir)
  exp_data <- read.csv(paste0("EXP_", cancername, ".csv"), row.names = 1)
  colnames(exp_data) <- chartr(old = '.', new = '_', colnames(exp_data))
  colnames(exp_data) <- chartr(old = '-', new = '_', colnames(exp_data))
  print(exp_data[1:5, 1:5])

  # ======================================================================
  
  setwd(cnv_dir)
  cnv_data <- read.csv(paste0("CNV_", cancername, ".csv"), row.names = 1)
  print(cnv_data[1:5, 1:5])
  colnames(cnv_data)[1:5]
  colnames(cnv_data) <- chartr(old = '.', new = '_', colnames(cnv_data))
  colnames(cnv_data) <- chartr(old = '-', new = '_', colnames(cnv_data))
  colnames(cnv_data)[1:5]
  
  metad <- data.frame(names(cnv_data))
  colnames(metad)[1] <- 'sample'
  metad$id <- substr(metad$sample, start = 1, stop = 12)
  metad$tape <- substr(metad$sample, start = 14, stop = 16)
  metad <- metad[order(metad$id), ]
  print(dim(metad))
  metad <- subset(metad, metad$tape == "01")
  
  print(dim(metad))
  metad <- metad %>% distinct(id, .keep_all = TRUE)
  print(dim(metad))
  
  cnv_data <- cnv_data[, which(colnames(cnv_data) %in% metad$sample)]
  colnames(cnv_data) <- substr(colnames(cnv_data), start = 1, stop = 12)
  cnv_data[1:5,1:5]
  
  min(cnv_data)
  
  # ======================================================================
  
  setwd(met_dir)
  met_data <- read.csv(paste0("Meth_", cancername, ".csv"), row.names = 1, na.strings = c(""))

  colnames(met_data) <- chartr(old = '.', new = '_', colnames(met_data))
  colnames(met_data) <- chartr(old = '-', new = '_', colnames(met_data))
  print(met_data[1:5, 1:5])

  
  metad <- data.frame(names(met_data))
  colnames(metad)[1] <- 'sample'
  metad$id <- substr(metad$sample, start = 1, stop = 12)
  metad$tape <- substr(metad$sample, start = 14, stop = 16)
  metad <- metad[order(metad$id), ]
  metad <- subset(metad, metad$tape == "01")
  
  print(dim(metad))
  metad <- metad %>% distinct(id, .keep_all = TRUE)
  print(dim(metad))
  
  met_data <- met_data[, which(colnames(met_data) %in% metad$sample)]
  colnames(met_data) <- substr(colnames(met_data), start = 1, stop = 12)
  met_data[1:5,1:5]
  

  # ======================================================================
  

  sele_data <- survivaldata[which(survivaldata$Type == cancername),]
  
  
  sample1 <- unique(sele_data$Sample)
  sample2 <- unique(colnames(exp_data))
  sample3 <- unique(colnames(cnv_data))
  sample4 <- unique(colnames(met_data))

  common_samples <- Reduce(intersect, list(sample1, sample2, sample3, sample4))
  
  
  exp_data <- exp_data[,which(colnames(exp_data) %in% common_samples)]
  cnv_data <- cnv_data[,which(colnames(cnv_data) %in% common_samples)]
  met_data <- met_data[,which(colnames(met_data) %in% common_samples)]

  sele_data <- sele_data[match(common_samples, sele_data$Sample), ]

  gene1 <- rownames(exp_data)
  gene2 <- rownames(cnv_data)
  gene3 <- rownames(met_data)

  common_genes <- Reduce(intersect, list(gene1, gene2, gene3))

  common_genes <- sort(common_genes)
  common_samples <- sort(common_samples)

  exp_data <- exp_data[common_genes, common_samples]
  cnv_data <- cnv_data[common_genes, common_samples]
  met_data <- met_data[common_genes, common_samples]
  
  
  setwd(save_dir)
  
  savefile <- paste0(cancername, "_EXP.csv")
  print(savefile)
  write.csv(exp_data, savefile)
  

  savefile <- paste0(cancername, "_CNV.csv")
  print(savefile)
  write.csv(cnv_data, savefile)

  savefile <- paste0(cancername, "_MET.csv")
  print(savefile)
  write.csv(met_data, savefile)

  savefile <- paste0(cancername, "_OS.csv")
  print(savefile)
  write.csv(sele_data, savefile)

  
  print(exp_data[1:5, 1:5])
  print(cnv_data[1:5, 1:5])
  print(met_data[1:5, 1:5])
  
  
  print('=============================================')
  
  
}
  



