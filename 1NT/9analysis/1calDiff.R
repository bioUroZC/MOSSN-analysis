#=============================================================

#=============================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(UpSetR)

data <- read.csv(
  paste0(PROJ_ROOT, "/1NT/2string/benchmark/MOSSN_uniform/merged_matrix.csv"),
  header = TRUE, row.names = 1, check.names = FALSE
)

cancers <- c("BLCA","BRCA","CRC","ESCA","HNSC",
             "KIRC","LIHC","LUAD","LUSC","PRAD","STAD")

setwd(paste0(PROJ_ROOT, "/1NT/9analysis/"))

cancertype <- read.csv(paste0(PROJ_ROOT, "/1NT/1data/TCGA/metadata.csv"),
                       header = TRUE, row.names = 1) %>%
  dplyr::rename(sample = Sample, cancer = Type) %>%
  dplyr::mutate(
    patient = base::substr(sample, 1, 12),
    group   = base::substr(sample, 14, 16),
    group   = dplyr::case_when(
      group == "01A" ~ "Cancer",
      group == "11A" ~ "Normal",
      TRUE           ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group))

res_list <- list()

#=============================================================

#=============================================================

for (cancername in cancers) {
  
  cancertype1 <- cancertype[cancertype$cancer %in% cancername, , drop = FALSE]
  
  common_samples <- intersect(colnames(data), cancertype1$sample)
  
  data_sub <- data[, common_samples, drop = FALSE]
  
  atlas_long <- data_sub %>%
    tibble::rownames_to_column("link") %>%
    tidyr::pivot_longer(
      cols = -link,
      names_to = "sample",
      values_to = "value"
    ) %>%
    dplyr::left_join(cancertype1, by = "sample") %>%   # IMPORTANT: use cancertype1
    dplyr::select(sample, link, group, patient, cancer, value)
  
  wide <- atlas_long %>%
    dplyr::select(patient, link, group, value) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      id_cols     = c(patient, link),
      names_from  = group,
      values_from = value
    ) %>%
    dplyr::filter(!base::is.na(Cancer) & !base::is.na(Normal)) %>%
    dplyr::mutate(diff = Cancer - Normal)
  
  atlas <- wide %>%
    dplyr::group_by(link) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      delta_median = stats::median(diff, na.rm = TRUE),
      delta_mean   = base::mean(diff, na.rm = TRUE),
      p_value = {
        x <- diff
        x <- x[!base::is.na(x)]
        if (base::length(x) == 0) NA_real_
        else if (base::length(base::unique(base::round(x, 6))) <= 1) 1
        else stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      fdr = stats::p.adjust(p_value, method = "BH")
    )
  
  ## tau and direction
  tau <- stats::quantile(base::abs(atlas$delta_median), 0.95, na.rm = TRUE)
  
  atlas <- atlas %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        !base::is.na(fdr) & fdr < 0.05 & delta_median >  tau ~ "gain",
        !base::is.na(fdr) & fdr < 0.05 & delta_median < -tau ~ "loss",
        TRUE ~ "no_change"
      ),
      tau_used = base::as.numeric(tau),
      cancer   = cancername
    ) %>%
    #dplyr::filter(direction %in% c("gain", "loss")) %>%
    base::as.data.frame()
  
  res_list[[cancername]] <- atlas
  
  base::message(
    cancername, ": total=", base::nrow(atlas),
    ", gain=", base::sum(atlas$direction == "gain", na.rm = TRUE),
    ", loss=", base::sum(atlas$direction == "loss", na.rm = TRUE)
  )
}

atlas_all <- dplyr::bind_rows(res_list)
head(atlas_all)

write.csv(atlas_all, file = "atlas_all.csv")
