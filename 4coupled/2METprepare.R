#=======================================================

#=======================================================

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#=======================================================
# Part 1: build probe -> gene promoter map (run once)
# Skip this block if promoter_probe_map.Rdata already exists
#=======================================================


rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


setwd(paste0(PROJ_ROOT, "/4coupled/METfiles"))

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- as.data.frame(anno)
probe_anno <- anno[, c("Name","chr","pos")]
probe_anno <- probe_anno[complete.cases(probe_anno), ]

probe_anno$chr <- ifelse(grepl("^chr", probe_anno$chr),
                         probe_anno$chr, paste0("chr", probe_anno$chr))

probe_gr <- GRanges(seqnames = probe_anno$chr,
                    ranges   = IRanges(start = probe_anno$pos,
                                       end = probe_anno$pos))
names(probe_gr) <- probe_anno$Name

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txs  <- transcripts(txdb, columns = c("tx_id","gene_id"))
prom_tx <- promoters(txs, upstream = 500, downstream = 50)

prom_tx$TXID   <- txs$tx_id
prom_tx$ENTREZ <- txs$gene_id

entrez_ids <- unique(unlist(prom_tx$ENTREZ))
entrez_ids <- entrez_ids[entrez_ids != ""]

sym <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = entrez_ids,
                             keytype = "ENTREZID",
                             columns = "SYMBOL")

prom_tx$SYMBOL <- sapply(prom_tx$ENTREZ, function(x) {
  if (length(x) == 0) return(NA)
  paste(unique(sym$SYMBOL[match(x, sym$ENTREZID)]), collapse = ";")
})

hits <- findOverlaps(probe_gr, prom_tx, ignore.strand = TRUE)
promoter_probe_map <- data.frame(
  probe = names(probe_gr)[queryHits(hits)],
  gene  = prom_tx$SYMBOL[subjectHits(hits)],
  stringsAsFactors = FALSE
)

promoter_probe_map <- unique(promoter_probe_map)
names(promoter_probe_map)[1] <- "Name"

anno1 <- merge(anno, promoter_probe_map, by="Name")
anno1 <- subset(anno1, select=c("UCSC_RefGene_Name", "gene"))

matches <- mapply(function(name, g) {
  if (is.na(g) || g == "") return(FALSE)
  grepl(paste0("\\b", g, "\\b"), name)
}, anno1$UCSC_RefGene_Name, anno1$gene)

sum(matches)
mean(matches)

save(promoter_probe_map, file = "promoter_probe_map.Rdata")

#=======================================================
# Part 2: process methylation data per cancer
#=======================================================

rm(list=ls())

setwd(paste0(PROJ_ROOT, "/4coupled/METfiles"))

load("promoter_probe_map.Rdata")

cancertypes <- c("ACC", "BLCA", "BRCA", "CESC", 
                 "CRC", "ESCA", "GBM", "HNSC", "KIRC",
                 "LGG", "LIHC", "LUAD", "LUSC", 
                 "PAAD", "PRAD", "SARC", "STAD")

for (i in 1:length(cancertypes)) {

  setwd(paste0(PROJ_ROOT, "/4coupled/METfiles"))

  infile <- paste0("TCGA.", cancertypes[i], ".sampleMap_HumanMethylation450.gz")
  print(cancertypes[i])

  methdata <- fread(infile)
  methdata <- as.data.frame(methdata)
  names(methdata)[1] <- "Name"

  na_rate <- rowMeans(is.na(methdata[,-1]))
  meth_filtered <- methdata[na_rate <= 0.3, ]

  mergeData <- merge(promoter_probe_map, meth_filtered, by="Name")
  mergeData$Name <- NULL

######use median to impute missing values for each probe
sample_cols <- setdiff(colnames(mergeData), c("Name", "gene"))
beta_mat <- as.matrix(mergeData[, sample_cols])

probe_medians <- apply(beta_mat, 1, median, na.rm = TRUE)
na_idx <- which(is.na(beta_mat), arr.ind = TRUE)
beta_mat[na_idx] <- probe_medians[na_idx[, 1]]
######

mergeData[, sample_cols] <- beta_mat

  exprSet <- aggregate(x = mergeData[, 2:ncol(mergeData)],
                       by = list(mergeData$gene),
                       FUN = mean)

  exprSet <- as.data.frame(exprSet)
  exprSet <- exprSet[-1, ]
  names(exprSet)[1] <- "ID"
  rownames(exprSet) <- exprSet$ID
  exprSet$ID <- NULL

  print(exprSet[1:5, 1:5])
  print(max(exprSet))
  print(min(exprSet))

  # invert beta values: high methylation -> low activity
  exprSet <- 1 - exprSet

  print(exprSet[1:5, 1:5])
  print(max(exprSet))
  print(min(exprSet))

  colnames(exprSet) <- chartr(old = "-", new = "_", colnames(exprSet))

  folder_for_save <- paste0(PROJ_ROOT, "/4coupled/METout/")
  dir.create(folder_for_save, showWarnings = FALSE)
  setwd(folder_for_save)

  file_for_save <- paste0("Meth_", cancertypes[i], ".csv")
  write.csv(exprSet, file = file_for_save)

}
