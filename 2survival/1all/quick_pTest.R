#=======================================================
#
#=======================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(officer)
library(flextable)
library(openxlsx)

#=======================================================
#
#=======================================================

cancertypes <- c("ACC", "BLCA", "BRCA", "CHOL",
                 "CRC", "GBM", "KIRC", "LGG",
                 "LIHC", "LUAD", "OV", "PAAD",
                 "PRAD",  "STAD")

method_levels <- c("MOSSN_uniform", "PPIXpress", "SSN",
                   "MOSSN_noCorr", "EdgeNoRWR", "RandomBackbone",
                   "RawExpr", "NodeRWR")

ref_method <- "MOSSN_uniform"

all_data <- list()
for (ct in cancertypes) {
  f <- paste0(paste0(PROJ_ROOT, "/2survival/"), ct, "/2quick/ml_dataset.csv")
  if (!file.exists(f)) { message("Skipping (not found): ", ct); next }
  df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  empty_cols <- which(colnames(df) == "")
  if (length(empty_cols) > 0)
    colnames(df)[empty_cols] <- paste0("RowID", seq_along(empty_cols))
  df$CancerType <- ct
  all_data[[ct]] <- df
}
data <- do.call(rbind, all_data)
data <- na.omit(data)

data$File[data$File == "PPPbi1"] <- "MOSSN_uniform"
data$File[data$File == "PPPbi2"] <- "EdgeNoRWR"
data$File[data$File == "PPPbi3"] <- "PPIXpress"
data$File[data$File == "PPPbi4"] <- "SSN"
data$File[data$File == "PPPbi5"] <- "RawExpr"
data$File[data$File == "PPPbi6"] <- "NodeRWR"
data$File[data$File == "PPPbi7"] <- "MOSSN_noCorr"
data$File[data$File == "PPPbi8"] <- "RandomBackbone"

data$File <- factor(data$File, levels = method_levels)

data <- dplyr::rename(data, AUC = Mean_tAUC)
data$AUC <- as.numeric(data$AUC)
data$C_index <- as.numeric(data$C_index)

methods <- levels(droplevels(data$File))

#=======================================================
# Paired Wilcoxon test: each method vs MOSSN_uniform (reference).
# key_cols is the pairing unit -- CancerType for the tumor-level
# table (one mean per tumor), CancerType+Dataset for the dataset-
# level table (no aggregation, every dataset is its own sample;
# note this treats datasets within a CancerType as independent --
# see quick_pTest_cluster.R for the cluster-aware alternatives).
#=======================================================

get_p_values <- function(df, key_cols, metric, colname) {
  map_df(methods, function(m) {
    df_m <- df %>%
      filter(File == m) %>%
      select(all_of(key_cols), value = all_of(metric))

    df_r <- df %>%
      filter(File == ref_method) %>%
      select(all_of(key_cols), ref_value = all_of(metric))

    df_pair <- inner_join(df_m, df_r, by = key_cols)

    vals_m <- df_pair$value
    vals_r <- df_pair$ref_value

    if (length(vals_m) == 0 || length(vals_r) == 0 ||
        all(is.na(vals_m)) || all(is.na(vals_r))) {
      p <- NA_real_
    } else if (all(vals_m == vals_r, na.rm = TRUE)) {
      p <- NA_real_
    } else {
      p <- suppressWarnings(wilcox.test(vals_m, vals_r, paired = TRUE)$p.value)
    }

    tibble(File = m, !!colname := p)
  })
}

format_p <- function(p) {
  if (is.na(p)) return("Ref")
  if (p < 0.01) return("<0.01")
  if (p < 0.05) return("<0.05")
  return(sprintf("%.2f", p))
}

build_summary <- function(df, key_cols) {
  p_cindex <- get_p_values(df, key_cols, "mean_C_index", "p_cindex")
  p_auc    <- get_p_values(df, key_cols, "mean_AUC",     "p_AUC")

  final_table <- df %>%
    group_by(File) %>%
    summarise(
      Mean_Cindex = mean(mean_C_index, na.rm = TRUE),
      SD_Cindex   = sd(mean_C_index,   na.rm = TRUE),
      Mean_AUC    = mean(mean_AUC,     na.rm = TRUE),
      SD_AUC      = sd(mean_AUC,       na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Cindex_mean_sd = paste0(round(Mean_Cindex, 3), "±", round(SD_Cindex, 3)),
      AUC_mean_sd    = paste0(round(Mean_AUC,    3), "±", round(SD_AUC,    3))
    ) %>%
    left_join(p_cindex, by = "File") %>%
    left_join(p_auc,    by = "File")

  final_table$p_cindex <- sapply(final_table$p_cindex, format_p)
  final_table$p_AUC    <- sapply(final_table$p_AUC,    format_p)

  final_table %>%
    mutate(File = factor(File, levels = method_levels)) %>%
    arrange(File) %>%
    dplyr::select(File, Cindex_mean_sd, p_cindex, AUC_mean_sd, p_AUC) %>%
    as.data.frame()
}

#=======================================================
# Tumor-level: one row per CancerType x File (mean over that
# tumor's datasets), n = number of cancer types
#=======================================================

df_tumor <- data %>%
  group_by(CancerType, File) %>%
  summarise(
    mean_C_index = mean(C_index, na.rm = TRUE),
    mean_AUC     = mean(AUC,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

table_tumor <- build_summary(df_tumor, "CancerType")
print(table_tumor)

#=======================================================
# Dataset-level: one row per CancerType+Dataset x File, no
# aggregation, n = number of datasets
#=======================================================

df_dataset <- data %>%
  select(CancerType, Dataset, File, mean_C_index = C_index, mean_AUC = AUC)

table_dataset <- build_summary(df_dataset, c("CancerType", "Dataset"))
print(table_dataset)

#=======================================================
#
#=======================================================

out_dir <- paste0(PROJ_ROOT, "/2survival/1all")
setwd(out_dir)

write.csv(table_tumor,   file = "quick_pTest_tumor.csv",   row.names = FALSE)
write.csv(table_dataset, file = "quick_pTest_dataset.csv", row.names = FALSE)

#=======================================================
# Merge both levels into a single table with a leading "Level"
# column ("Tumor" / "Dataset"), written to one .docx and one .xlsx
#=======================================================

table_combined <- bind_rows(
  table_tumor   %>% mutate(Level = "Tumor",   .before = 1),
  table_dataset %>% mutate(Level = "Dataset", .before = 1)
)

ft <- flextable(table_combined)
ft <- autofit(ft)
ft <- theme_vanilla(ft)
ft <- align(ft, align = "center", part = "all")
ft <- fontsize(ft, size = 10, part = "all")
ft <- bold(ft, part = "header")
ft <- merge_v(ft, j = "Level")

doc <- read_docx()
doc <- body_add_par(doc, "Survival model comparison (2quick, vs MOSSN_uniform)", style = "heading 1")
doc <- body_add_par(doc, "")
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "quick_pTest.docx")

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", table_combined)
saveWorkbook(wb, "quick_pTest.xlsx", overwrite = TRUE)

cat("quick_pTest_tumor.csv / quick_pTest_dataset.csv, quick_pTest.docx, quick_pTest.xlsx written to", out_dir, "\n")
