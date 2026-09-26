rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(officer)
library(flextable)

cancertypes <- c("ACC", "BLCA", "BRCA", "CHOL",
                 "CRC", "GBM", "KIRC", "LGG",
                 "LIHC", "LUAD", "OV", "PAAD",
                 "PRAD",  "STAD")

method_levels <- c("MOSSN", "PPIXpress", "Patkar",
                   "noCorr", "noRWR", "RandomBackbone",
                   "RawExpr", "NodeRWR")

ref_method <- "MOSSN"

all_data <- list()
for (ct in cancertypes) {
  f <- paste0("/proj/c.zihao/work1/03survival/", ct, "/03sur/ml_dataset.csv")
  if (!file.exists(f)) {
    message("Skipping (not found): ", ct)
    next
  }
  df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  empty_cols <- which(colnames(df) == "")
  if (length(empty_cols) > 0)
    colnames(df)[empty_cols] <- paste0("RowID", seq_along(empty_cols))
  df$CancerType <- ct
  all_data[[ct]] <- df
}
data <- do.call(rbind, all_data)
data <- na.omit(data)

data$File <- factor(data$File, levels = method_levels)

data <- dplyr::rename(data, AUC = Mean_tAUC)
data$AUC <- as.numeric(data$AUC)
data$C_index <- as.numeric(data$C_index)

methods <- levels(droplevels(data$File))

get_p_values <- function(df, key_cols, metric, colname) {
  map_df(methods, function(m) {
    df_m <- df %>%
      dplyr::filter(File == m) %>%
      dplyr::select(all_of(key_cols), value = all_of(metric))

    df_r <- df %>%
      dplyr::filter(File == ref_method) %>%
      dplyr::select(all_of(key_cols), ref_value = all_of(metric))

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
    dplyr::select(
      Method = File,
      `C-index (mean±sd)` = Cindex_mean_sd,
      `p-value (C-index)` = p_cindex,
      `timeAUC (mean±sd)` = AUC_mean_sd,
      `p-value (timeAUC)` = p_AUC
    ) %>%
    as.data.frame()
}

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

df_dataset <- data %>%
  dplyr::select(CancerType, Dataset, File, mean_C_index = C_index, mean_AUC = AUC)

table_dataset <- build_summary(df_dataset, c("CancerType", "Dataset"))
print(table_dataset)

out_dir <- "/proj/c.zihao/work1/03survival/1all"
setwd(out_dir)

table_combined <- bind_rows(
  table_tumor   %>% mutate(Level = "Tumor",   .before = 1),
  table_dataset %>% mutate(Level = "Dataset", .before = 1)
)

write.csv(table_combined, file = "Table 3.csv", row.names = FALSE)

ft <- flextable(table_combined)
ft <- autofit(ft)
ft <- theme_vanilla(ft)
ft <- align(ft, align = "center", part = "all")
ft <- fontsize(ft, size = 10, part = "all")
ft <- bold(ft, part = "header")
ft <- merge_v(ft, j = "Level")

doc <- read_docx()
doc <- body_add_par(
  doc,
  paste0(
    "Table 3. Survival prediction performance of MOSSN compared with external ssPIN methods and MOSSN ablation variants. ",
    "Performance was evaluated under leave-one-dataset-out cross-validation using LASSO-Cox regression, reported as the C-index and time-dependent AUC (timeAUC; mean ± sd) at the tumor level (averaged within each cancer type) and the dataset level (each dataset treated independently). ",
    "P-values are from paired Wilcoxon signed-rank tests against MOSSN (Ref)."
  ),
  style = "Normal"
)
doc <- body_add_par(doc, "")
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "Table 3.docx")

cat("Table 3.csv / Table 3.docx written to", out_dir, "\n")
