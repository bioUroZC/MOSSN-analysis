rm(list=ls())

setwd("/proj/c.zihao/work1/01data/01count")

library(readxl)
library(dplyr)

df <- read_excel("dataset.xlsx", sheet = "Sheet1")

n_datasets <- nrow(df)
n_samples <- sum(df[["Sample Number"]], na.rm = TRUE)

cat("Total datasets:", n_datasets, "\n")
cat("Total samples:", n_samples, "\n\n")

by_cancer <- df %>%
  group_by(`Cancer Type`) %>%
  summarise(
    n_datasets = n(),
    n_samples = sum(`Sample Number`, na.rm = TRUE),
    .groups = "drop"
  )

cat("By cancer type:\n")
print(by_cancer, n = Inf)

tech <- case_when(
  grepl("HiSeq|NextSeq|NovaSeq|Genome Analyzer|X Ten|RNA-Seq", df$Platform, ignore.case = TRUE) ~ "RNA-seq",
  grepl("Array|beadchip|Microarray|DASL", df$Platform, ignore.case = TRUE) ~ "Microarray",
  TRUE ~ "Unclassified"
)

cat("\nMicroarray datasets:", sum(tech == "Microarray"), "\n")
cat("RNA-seq datasets:", sum(tech == "RNA-seq"), "\n")
cat("Unclassified datasets:", sum(tech == "Unclassified"), "\n")
print(unique(df$Platform[tech == "Unclassified"]))
