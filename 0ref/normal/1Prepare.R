#=======================================================

#=======================================================


rm(list=ls())
options(stringsAsFactors = F)
library(biomaRt)
library(dplyr)

#=======================================================

#=======================================================

rm(list=ls())

setwd("/proj/c.zihao/work1/0ref/GTEx/")
rdseed <- 8

load('GTEx_all.Rdata')
GTEx[1:5,1:5]
colnames(GTEx) <-  gsub('[.]','-',colnames(GTEx)) 


dim(GTEx)
GTEx <- GTEx[!duplicated(GTEx$Description), ]
rownames(GTEx) <- GTEx$Description
dim(GTEx)
GTEx$Name <- NULL
GTEx$Description  <- NULL
GTEx[1:5,1:5]
# GTEx <- log2(GTEx+1)
print(max(GTEx))
print(min(GTEx))

#=======================================================

# Connect to Ensembl (use GRCh38/hg38; adjust if using GRCh37)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Fetch gene names and their biotype
gene_info <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  mart = mart
)

head(gene_info)

# Filter the annotation table
protein_coding_genes <- gene_info %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::pull(external_gene_name)


# Keep only rows (genes) that are protein coding
GTEx <- GTEx[rownames(GTEx) %in% protein_coding_genes, ]

#=======================================================

pheno <- read.table('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep = ' ', header = T)

pheno <-  subset(pheno,  select = c( "SAMPID", "SMTS" ))
dt <-  as.data.frame( table(pheno$SMTS) )

dt

seleTissu <- c("Adrenal Gland",
               'Bladder',	
               'Brain',	
               'Breast',	
               'Cervix Uteri',	
               'Colon',	
               'Esophagus',	
               'Kidney',	
               'Liver',	
               'Lung',
               "Muscle",
               'Ovary',	
               'Pancreas',	
               'Prostate',
               'Skin',
               'Stomach',
               "Salivary Gland")	

pheno <- pheno[which(pheno$SMTS %in% seleTissu),]

head(pheno)

dt <-  as.data.frame(table(pheno$SMTS) )

dt

#=======================================================

# Initialize list to store sample-wise dataframes
sample_df_list <- list()

# Unique organs
organs <- unique(pheno$SMTS)

set.seed(rdseed)

for (organ in organs) {
  
  print(organ)
        
  sample_meta <- pheno %>%  dplyr::filter(SMTS == organ)
  organ_expr <- GTEx[, colnames(GTEx) %in% sample_meta$SAMPID]
  
  # Sample 10 columns if available
  if (ncol(organ_expr) >= 10) {
    sampled_cols <- sample(ncol(organ_expr), 10)
    organ_expr <- organ_expr[, sampled_cols]
  }
  
  # Transpose expression data: samples become rows
  organ_expr_t <- as.data.frame(t(organ_expr))
  
  # Add sample and organ info
  organ_expr_t <- organ_expr_t %>%
     dplyr::mutate(
      sample_id = rownames(organ_expr_t),
      organ = organ
    ) %>%
     dplyr::relocate(sample_id, organ)  # Move metadata columns to front
  
  sample_df_list[[organ]] <- organ_expr_t
}

# Combine all organs into one dataframe
combined_expr_df <- dplyr::bind_rows(sample_df_list)
combined_expr_df[1:5,1:5]
rownames(combined_expr_df) <- combined_expr_df$sample_id
combined_expr_df$sample_id <- NULL

save_path <- '/proj/c.zihao/work1/0ref/normal/'
print(save_path)

setwd(save_path)
write.csv(combined_expr_df, file = "combined_expr_df.csv")
print("All finished")
