table_file <- "/proj/c.zihao/work1/03survival/0sup/Supplementary_Table_2.csv"
output_file <- "/proj/c.zihao/work1/03survival/0sup/expression_modality_summary.csv"

dat <- read.csv(table_file, check.names = FALSE)

is_rna_seq <- grepl("RNA-seq|HiSeq|NovaSeq|NextSeq",
                     dat$Platform,
                     ignore.case = TRUE)

dat$Modality <- "Microarray"
dat$Modality[is_rna_seq] <- "RNA-seq"

microarray <- dat[dat$Modality == "Microarray", ]
rna_seq <- dat[dat$Modality == "RNA-seq", ]

result <- data.frame(
  Modality = c("Microarray", "RNA-seq", "Total"),
  `Dataset Number` = c(nrow(microarray), nrow(rna_seq), nrow(dat)),
  `Sample Number` = c(
    sum(microarray$`Sample Number`),
    sum(rna_seq$`Sample Number`),
    sum(dat$`Sample Number`)
  ),
  check.names = FALSE
)

write.csv(result, output_file, row.names = FALSE)
print(result)
