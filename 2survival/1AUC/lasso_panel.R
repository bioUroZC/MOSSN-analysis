rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)

cancertypes <- c("ACC", "BLCA", "BRCA", "CHOL",
                 "CRC",   "GBM", "KIRC",
                 "LGG",  "LIHC", "LUAD", "OV",
                 "PAAD", "PRAD",  "STAD")

method_levels <- c("MOSSN_noPrior", "MOSSN_uniform", "Patkar", "PPIXpress",
                   "Proteinarium", "SSN", "SWEET", "LIONESS")

load_data <- function(suffix) {
  all_data <- list()
  for (ct in cancertypes) {
    f <- paste0("/proj/c.zihao/work1/2survival/", ct, "/2sur/", suffix, ".csv")
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
  data$File[data$File == "PPPbi1"] <- "MOSSN_noPrior"
  data$File[data$File == "PPPbi2"] <- "MOSSN_uniform"
  data$File[data$File == "PPPbi3"] <- "Patkar"
  data$File[data$File == "PPPbi4"] <- "PPIXpress"
  data$File[data$File == "PPPbi5"] <- "Proteinarium"
  data$File[data$File == "PPPbi6"] <- "SSN"
  data$File[data$File == "PPPbi7"] <- "SWEET"
  data$File[data$File == "PPPbi8"] <- "LIONESS"
  data$File <- factor(data$File, levels = method_levels)
  data
}

summarise_data <- function(data) {
  dt <- data %>%
    group_by(CancerType, File) %>%
    summarise(
      mean_C_index = mean(C_index,   na.rm = TRUE),
      mean_tAUC    = mean(Mean_tAUC, na.rm = TRUE),
      .groups = "drop"
    )
  dt$CancerType <- factor(dt$CancerType, levels = cancertypes)
  dt <- na.omit(dt)

  overall_means <- dt %>%
    group_by(File) %>%
    summarise(
      overall_mean_C_index = mean(mean_C_index, na.rm = TRUE),
      overall_mean_tAUC    = mean(mean_tAUC,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(overall_mean_tAUC))

  list(dt = dt, overall_means = overall_means)
}

out_dir <- "/proj/c.zihao/work1/2survival/1AUC"

dataset_res <- summarise_data(load_data("lasso_dataset"))
rep_res     <- summarise_data(load_data("lasso_rep"))

write.csv(dataset_res$dt,            file.path(out_dir, "lasso_dataset_tumor_mean.csv"),    row.names = FALSE)
write.csv(dataset_res$overall_means, file.path(out_dir, "lasso_dataset_overall_means.csv"), row.names = FALSE)
write.csv(rep_res$dt,                file.path(out_dir, "lasso_rep_tumor_mean.csv"),         row.names = FALSE)
write.csv(rep_res$overall_means,     file.path(out_dir, "lasso_rep_overall_means.csv"),      row.names = FALSE)

print(dataset_res$overall_means)
print(rep_res$overall_means)

#=======================================================

custom_colors <- c(
  MOSSN_noPrior = "#67000d",
  MOSSN_uniform = "#a50f15",
  SSN           = "#4c78a8",
  SWEET         = "#72b7b2",
  Patkar        = "#c9a227",
  Proteinarium  = "#9d755d",
  PPIXpress     = "#b279a2",
  LIONESS       = "#54a24b"
)

make_bar <- function(dt, y_var, y_label, ylim_top = NULL) {
  p <- ggplot(dt, aes(x = CancerType, y = .data[[y_var]], fill = File)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Cancer Type", y = y_label, fill = "Method") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_manual(values = custom_colors)
  if (!is.null(ylim_top)) p <- p + ylim(0.0, ylim_top)
  p
}

p1 <- make_bar(dataset_res$dt, "mean_tAUC",    "Mean tAUC (Dataset)",    0.9)
p2 <- make_bar(rep_res$dt,     "mean_tAUC",    "Mean tAUC (Rep)",        0.9)
p3 <- make_bar(dataset_res$dt, "mean_C_index", "Mean C-index (Dataset)")
p4 <- make_bar(rep_res$dt,     "mean_C_index", "Mean C-index (Rep)")

panel <- (p1 | p3) / (p2 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file.path(out_dir, "lasso_panel.pdf"), height = 10, width = 16)
print(panel)
dev.off()
