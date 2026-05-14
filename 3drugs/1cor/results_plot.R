rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)

cancertypes <- c("CisplatinSTAD", "FluoroBRCA", "imMelanoma", "Paclitaxel")

method_levels <- c("MOSSN", "MOSSN (uniform)", "Patkar", "PPIXpress",
                   "Proteinarium", "SSN", "SWEET", "LIONESS")


custom_colors <- c(
  "MOSSN"           = "#67000d",
  "MOSSN (uniform)" = "#a50f15",
  "SSN"             = "#4c78a8",
  "SWEET"           = "#72b7b2",
  "Patkar"          = "#c9a227",
  "Proteinarium"    = "#9d755d",
  "PPIXpress"       = "#b279a2",
  "LIONESS"         = "#54a24b"
)

load_data <- function(file_suffix) {
  all_data <- list()
  for (ct in cancertypes) {
    f <- paste0("/proj/c.zihao/work1/3drugs/", ct, "/2response/", file_suffix)
    if (!file.exists(f)) { message("Skipping (not found): ", ct); next }
    df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    empty_cols <- which(colnames(df) == "")
    if (length(empty_cols) > 0)
      colnames(df)[empty_cols] <- paste0("RowID", seq_along(empty_cols))
    df$CancerType <- ct
    all_data[[ct]] <- df
  }
  data <- do.call(rbind, all_data)
  data$File[data$File == "PPPbi1"] <- "MOSSN"
  data$File[data$File == "PPPbi2"] <- "MOSSN (uniform)"
  data$File[data$File == "PPPbi3"] <- "Patkar"
  data$File[data$File == "PPPbi4"] <- "PPIXpress"
  data$File[data$File == "PPPbi5"] <- "Proteinarium"
  data$File[data$File == "PPPbi6"] <- "SSN"
  data$File[data$File == "PPPbi7"] <- "SWEET"
  data$File[data$File == "PPPbi8"] <- "LIONESS"
  data$File <- factor(data$File, levels = method_levels)
  data$CancerType <- factor(data$CancerType, levels = cancertypes)
  data
}

summarise_data <- function(data) {
  dt <- data %>%
    group_by(CancerType, File) %>%
    summarise(mean_Spearman = mean(Spearman, na.rm = TRUE), .groups = "drop")
  overall <- dt %>%
    group_by(File) %>%
    summarise(overall_mean_Spearman = mean(mean_Spearman, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(overall_mean_Spearman))
  list(dt = dt, overall = overall)
}

make_plot <- function(dt, title) {
  ggplot(dt, aes(x = CancerType, y = mean_Spearman, fill = File)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Cancer Type", y = "Mean Spearman", fill = "Method",
         title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_manual(values = custom_colors)
}

# Dataset
res_dataset <- summarise_data(load_data("results_dataset.csv"))
cat("=== Dataset ===\n"); print(res_dataset$overall)

# Rep
res_rep <- summarise_data(load_data("results_rep.csv"))
cat("=== Rep ===\n"); print(res_rep$overall)

out_dir <- "/proj/c.zihao/work1/3drugs/1cor"
write.csv(res_dataset$overall,
          file.path(out_dir, "results_dataset_overall_means.csv"),
          row.names = FALSE)
write.csv(res_rep$overall,
          file.path(out_dir, "results_rep_overall_means.csv"),
          row.names = FALSE)

p1 <- make_plot(res_dataset$dt, "Dataset")
p2 <- make_plot(res_rep$dt, "Rep")

combined <- p1 + p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file.path(out_dir, "results_combined_Spearman.pdf"), height = 6, width = 14)
print(combined)
dev.off()
