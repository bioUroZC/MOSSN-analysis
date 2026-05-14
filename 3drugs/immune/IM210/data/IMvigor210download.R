# ===================================================

# ===================================================

rm(list = ls())

library(DESeq2)
library(DESeq)

# install.packages("IMvigor210CoreBiologies_1.0.1.tar.gz", 
#                  repos=NULL)

library("IMvigor210CoreBiologies")

data(cds)

head(counts(cds))[1:5,1:5]
head(fData(cds))[1:5,1:5]
head(pData(cds))[1:5,1:5]


expr <- data.frame(counts(cds))
anno <- data.frame(fData(cds))

anno <- subset(anno, select=c('entrez_id', 'symbol'))

expr$entrez_id <- rownames(expr)

data <- merge(anno, expr, by="entrez_id")

data[1:5,1:5]

data$entrez_id <- NULL

data[1:5,1:5]

exprSet <- aggregate(x = data[,2:ncol(data)],
                     by = list(data$symbol),
                     FUN = max)
exprSet[1:5,1:5]

exprSet <- exprSet[-1,]

exprSet[1:5,1:5]

names(exprSet)[1] <- 'gene'

exprSet[1:5,1:5]

write.csv(exprSet, file = "IMvigor210_exprSet.csv", row.names = F)


pheno <- data.frame(pData(cds))

head(pheno)

write.csv(pheno, "IMvigor210_FollowUp.csv")



