#===================================================

#===================================================

rm(list = ls())
library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)  

#===================================================

#===================================================

dieasename <- "Fluorouracil"
projecname <- "GSE20194"

base_path <- paste0("/proj/c.zihao/work1/3drugs/",dieasename, "/", projecname)
base_path 

sifout_file <- paste0(base_path, "/PPIX/human_ppin.sif.gz")
sifout_file

setwd('/proj/c.zihao/work1/1survival')
links <- read.csv("links.csv", row.names = 1)
head(links)

links <- links %>%
  dplyr::select(protein1, protein2, score)

names(links) <- c("Protein1", "Protein2", "weight")
head(links)

genes <- unique(union(links$Protein1, links$Protein2))

converted_genes <- bitr(genes,
                        fromType = "SYMBOL",
                        toType = "UNIPROT",
                        OrgDb = org.Hs.eg.db)

print(head(converted_genes))
head(links)

links_uniprot <- merge(links, converted_genes, 
                       by.x = "Protein1", 
                       by.y = "SYMBOL", all.x = TRUE)

names(links_uniprot)[ncol(links_uniprot)] <- "UniProt1"
head(links_uniprot)

links_uniprot <- merge(links_uniprot, converted_genes, by.x = "Protein2",
                       by.y = "SYMBOL", all.x = TRUE)

head(links_uniprot)
names(links_uniprot)

names(links_uniprot)[ncol(links_uniprot)] <- "UniProt2"
head(links_uniprot)

links_uniprot <- links_uniprot[, c("UniProt1", "UniProt2", "weight")]
head(links_uniprot)

names(links_uniprot) <- c('Protein1', 'Protein2', 'weight' )

links_uniprot$weight <- 1

write.table(links_uniprot, gzfile(sifout_file), 
            row.names = FALSE, col.names = T, sep = "\t", quote = FALSE)

cat("File saved as", sifout_file)

#===================================================

#===================================================

wddir <- paste0(base_path, '/data')
wddir 
setwd(wddir)

dtaa <- read.csv("exprSet_filtered.csv", header = T, row.names = 1)
dtaa[1:5,1:5] 
min(dtaa)
max(dtaa)

thresoldNumn <- median(as.numeric(unlist(dtaa)))
thresoldNumn
thresoldNumn <- round(thresoldNumn)
thresoldNumn

str(dtaa)

dtaa$SYMBOL <- rownames(dtaa)

gene_list <- rownames(dtaa)

converted_ids <- bitr(
  gene_list, 
  fromType = "SYMBOL",
  toType = c("ENSEMBL"),
  OrgDb = org.Hs.eg.db)

datase <- merge(converted_ids, dtaa, by="SYMBOL")

datase$SYMBOL <- NULL

print(head(datase)[1:5,1:5])


save_path <- file.path(base_path, "/PPIX/files/")
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)


for (p in 2:dim(datase)[2]) {
  
  selesamples <- colnames(datase)[p]
  print(selesamples)
  datase1 <- datase[,c('ENSEMBL',selesamples)]
  datase1 <- as.data.frame(datase1)
  names(datase1)[2] <- "normalized_count"
  
  file_path <- paste0(save_path, selesamples, ".gz")
  print(file_path)
  
  write.table(datase1, gzfile(file_path), row.names = FALSE, 
              col.names = FALSE, sep = "\t", quote = FALSE)
  
}


PPINout_path <- file.path(base_path, "/PPIX/out/")
dir.create(PPINout_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(PPINout_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)


#===================================================


#===================================================

data_dir <- paste0(base_path, "/PPIX/files")
data_dir

jar_path <- paste0(base_path, "/PPIX/PPIXpress123")
jar_path

sif_file <- paste0(base_path, "/PPIX/human_ppin.sif.gz")
sif_file 

# Get all .gz files
gz_files <- list.files(
  path = data_dir,
  pattern = "\\.gz$",
  full.names = TRUE
)

gz_files


# Count files
num_gz_files <- length(gz_files)
cat("Number of .gz files found:", num_gz_files, "\n")


output_dir <- paste0(base_path, "/PPIX/out")
output_dir  

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

thresoldNumn <- gsub(" ", "", paste("-t=", thresoldNumn))
thresoldNumn


# Build Java command
java_cmd <- paste(
  "java -jar PPIXpress.jar", thresoldNumn, sif_file, output_dir,
  paste(gz_files, collapse = " \\\n")
)

# Slurm script content
slurm_script <- c(
  "#!/bin/bash",
  "#SBATCH -J 1PPP",
  "#SBATCH -N 1",
  "#SBATCH --ntasks-per-node=1",
  "#SBATCH --mem=40G",
  "#SBATCH -o 1PPP.out",
  "#SBATCH -e 1PPP.err",
  "",
  paste0("cd ", jar_path),
  java_cmd
)

# Write script to file
script_path <- file.path(jar_path, paste0("run_PPP.sh"))
writeLines(slurm_script, script_path)
Sys.chmod(script_path, "0755")

cat("Generated Slurm script:", script_path, "\n")
thresoldNumn
