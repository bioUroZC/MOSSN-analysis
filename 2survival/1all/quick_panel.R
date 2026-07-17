rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(ggplot2)
library(patchwork)

cancertypes <- c("ACC", "BLCA", "BRCA", "CHOL",
                 "CRC", "GBM", "KIRC", "LGG",
                 "LIHC", "LUAD", "OV", "PAAD",
                 "PRAD", "STAD")

method_levels <- c("MOSSN_uniform", "PPIXpress", "SSN",
                   "MOSSN_noCorr", "EdgeNoRWR", "RandomBackbone",
                   "RawExpr", "NodeRWR")

load_data <- function() {
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

out_dir <- paste0(PROJ_ROOT, "/2survival/1all")

raw_data <- load_data()
quick_res <- summarise_data(raw_data)

write.csv(quick_res$dt,            file.path(out_dir, "quick_tumor_mean.csv"),    row.names = FALSE)
write.csv(quick_res$overall_means, file.path(out_dir, "quick_overall_means.csv"), row.names = FALSE)

print(quick_res$overall_means)

# Significance testing (MOSSN_uniform vs each method, paired Wilcoxon
# across CancerType) lives in quick_pTest.R, following the same
# convention as 1AUCsum/3setTest.R.

#=======================================================

custom_colors <- c(
  MOSSN_uniform  = "#67000d",
  EdgeNoRWR      = "#4c78a8",
  PPIXpress      = "#b279a2",
  SSN            = "#f58518",
  RawExpr        = "#72b7b2",
  NodeRWR        = "#54a24b",
  MOSSN_noCorr   = "#9d755d",
  RandomBackbone = "#bab0ac"
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

p1 <- make_bar(quick_res$dt, "mean_tAUC",    "Mean tAUC",    0.9)
p2 <- make_bar(quick_res$dt, "mean_C_index", "Mean C-index", 0.9)

panel <- (p1 / p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf(file.path(out_dir, "quick_panel.pdf"), height = 8, width = 8)
print(panel)
dev.off()

#=======================================================
# Dataset-level boxplot: each point is one dataset (pooled across
# all cancer types), x-axis is the 8 methods.
#=======================================================

lighten_color <- function(color, amount = 0.55) {
  rgb_val <- grDevices::col2rgb(color) / 255
  rgb_val <- rgb_val + (1 - rgb_val) * amount
  grDevices::rgb(rgb_val[1, ], rgb_val[2, ], rgb_val[3, ])
}

point_colors <- setNames(lighten_color(custom_colors), names(custom_colors))

make_box <- function(df, y_var, y_label, ylim_range = NULL) {
  p <- ggplot(df, aes(x = File, y = .data[[y_var]], fill = File)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = File), width = 0.15, size = 1.2, alpha = 0.8) +
    labs(x = "Method", y = y_label, fill = "Method") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = point_colors, guide = "none")
  if (!is.null(ylim_range)) p <- p + coord_cartesian(ylim = ylim_range)
  p
}

b1 <- make_box(raw_data, "Mean_tAUC", "Mean tAUC",   ylim_range = c(0.3, 0.8))
b2 <- make_box(raw_data, "C_index",   "Mean C-index", ylim_range = c(0.3, 0.8))

box_panel <- (b1 | b2) &
  theme(legend.position = "none")

pdf(file.path(out_dir, "quick_boxplot.pdf"), height = 5, width = 14)
print(box_panel)
dev.off()
