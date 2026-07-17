#=======================================================

#=======================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(data.table)

setwd(paste0(PROJ_ROOT, '/1NT/2CPTAC/eLUAD'))

#=======================================================

#=======================================================

dat <- fread(paste0(PROJ_ROOT, '/1NT/2CPTAC/eLUAD/GSE229705_counts-raw.csv.gz'))
dat <- as.data.frame(dat)
print(dat[1:5, 1:5])

all_samples <- colnames(dat)[-1]
tumor_samples  <- all_samples[grepl("-T$", all_samples)]
normal_samples <- all_samples[grepl("-N$", all_samples)]

tumor_patients  <- sub("-T$", "", tumor_samples)
normal_patients <- sub("-N$", "", normal_samples)
paired_patients <- intersect(tumor_patients, normal_patients)

cat("Tumor samples:", length(tumor_samples),
    "\nNormal samples:", length(normal_samples),
    "\nPaired patients:", length(paired_patients), "\n")

tumor_paired  <- dat[, c("gene", paste0(paired_patients, "-T")), drop = FALSE]
normal_paired <- dat[, c("gene", paste0(paired_patients, "-N")), drop = FALSE]

colnames(tumor_paired)[-1]  <- paste0(paired_patients, "_Tumor")
colnames(normal_paired)[-1] <- paste0(paired_patients, "_Normal")

mergedData <- merge(tumor_paired, normal_paired, by = "gene")

colnames(mergedData) <- gsub("-", "_", colnames(mergedData))
mergedData[is.na(mergedData)] <- 0

dim(mergedData)
print(mergedData[1:5, 1:4])

rownames(mergedData) <- mergedData$gene
mergedData$gene <- NULL

link_df <- fread(paste0(PROJ_ROOT, "/1NT/1data/string/links.csv"))
link_genes <- unique(c(as.character(link_df$protein1), as.character(link_df$protein2)))
keep_genes <- intersect(rownames(mergedData), link_genes)
mergedData <- mergedData[keep_genes, , drop = FALSE]

cat("Genes retained after STRING link filtering:", nrow(mergedData), "\n")

print(min(mergedData))
print(max(mergedData))

mergedData <- log2(mergedData + 1)


print(min(mergedData))
print(max(mergedData))

write.csv(mergedData, file = "eLUAD_exprSet_filtered.csv", row.names = TRUE)
print(dim(mergedData))

print(mergedData[1:5,1:5])
print(colnames(mergedData))
