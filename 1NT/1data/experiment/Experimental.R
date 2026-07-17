# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)
library(igraph)
library(poweRlaw)

# ===================================================

# ===================================================

setwd(paste0(PROJ_ROOT, "/2survival"))
protein_links <- fread("9606.protein.links.detailed.v12.0.txt", header = TRUE)
protein_links <- as.data.frame(protein_links)
protein_links <- subset(protein_links, select=c('protein1', 'protein2', 'experimental'))
str(protein_links)

min(protein_links$experimental)
max(protein_links$experimental)

protein_links <- protein_links %>%
  dplyr::filter(experimental > 900)
print(head(protein_links))


# ===================================================

# ===================================================

geneinfor <- fread("9606.protein.aliases.v12.0.txt")
geneinfor  <- as.data.frame(geneinfor)
table(geneinfor$source)
hugo_rows <- geneinfor[grep("HUGO", geneinfor$source), ]
print(head(hugo_rows))
names(hugo_rows)[1] <- "protein_id"


id_to_alias <- setNames(hugo_rows$alias, hugo_rows$protein_id)
protein_links$protein1_alias <- id_to_alias[protein_links$protein1]
protein_links$protein2_alias <- id_to_alias[protein_links$protein2]


print(head(protein_links))

protein_links$protein1 <- NULL
protein_links$protein2 <- NULL


head(protein_links)

# ===================================================

# ===================================================


# Ensure undirected uniqueness by alphabetically sorting alias pairs
protein_links <- protein_links %>%
  mutate(
    node1 = pmin(protein1_alias, protein2_alias),
    node2 = pmax(protein1_alias, protein2_alias)
  ) %>%
  dplyr::select(experimental, node1, node2) 

head(protein_links)

names(protein_links) <- c("score", "protein1", "protein2")

head(protein_links)

dim(protein_links)

protein_links <- protein_links %>%
  filter(protein1 != protein2)

dim(protein_links)

protein_links$link <- paste0(protein_links$protein1, '_', protein_links$protein2)

print(dim(protein_links))

protein_links <- protein_links %>%
  dplyr::distinct(link, .keep_all = T)

protein_links$link <- NULL

str(protein_links)

print(dim(protein_links))

protein_links$score <- protein_links$score / 1000

setwd(paste0(PROJ_ROOT, "/1NT/1data/experiment"))

write.csv(protein_links, file = 'experimental_link.csv')


