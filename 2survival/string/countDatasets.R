library(dplyr)

df <- read.csv("surDatasets.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(df) <- trimws(colnames(df))
df$`Sample Number` <- as.numeric(trimws(df$`Sample Number`))

# Sequencing: HiSeq [number], NovaSeq, NextSeq [number], Genome Analyzer, HiSeq X
is_seq <- grepl("HiSeq [0-9]|NovaSeq|NextSeq [0-9]|Genome Analyzer|HiSeq X", df$Platform)

df$DataType <- ifelse(is_seq, "Sequencing", "Microarray")

summary_stats <- df %>%
  group_by(DataType) %>%
  summarise(
    n_datasets = n(),
    n_samples  = sum(`Sample Number`, na.rm = TRUE),
    .groups = "drop"
  )

total <- tibble(
  DataType   = "Total",
  n_datasets = nrow(df),
  n_samples  = sum(df$`Sample Number`, na.rm = TRUE)
)

result <- bind_rows(summary_stats, total)
print(result, n = Inf)
