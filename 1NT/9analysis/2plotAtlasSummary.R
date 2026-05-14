rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
setwd(base_dir)

input_file <- file.path(base_dir, "atlas_all.csv")
out_gain_pdf <- file.path(base_dir, "gain_loss.pdf")
out_pair_pdf <- file.path(base_dir, "pairwise_heatmap.pdf")

cancers <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
             "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

atlas_all <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("" %in% names(atlas_all)) {
  atlas_all <- atlas_all[, names(atlas_all) != "", drop = FALSE]
}

plot_gain_loss_by_cancer <- function(atlas_df, cancer_levels, out_pdf) {
  plot_df <- atlas_df %>%
    dplyr::filter(direction %in% c("gain", "loss"), cancer %in% cancer_levels) %>%
    dplyr::count(cancer, direction, name = "n_links") %>%
    tidyr::complete(
      cancer = cancer_levels,
      direction = c("gain", "loss"),
      fill = list(n_links = 0L)
    ) %>%
    dplyr::mutate(
      cancer = factor(cancer, levels = cancer_levels),
      direction = factor(direction, levels = c("gain", "loss")),
      y = dplyr::if_else(direction == "loss", -n_links, n_links)
    )
  
  max_n <- max(plot_df$n_links, na.rm = TRUE)
  offset <- max_n * 0.02
  
  p <- ggplot(plot_df, aes(x = cancer, y = y, fill = direction)) +
    geom_col(width = 0.72) +
    geom_text(
      aes(
        label = n_links,
        y = ifelse(direction == "gain", y + offset, y - offset)
      ),
      size = 3.2
    ) +
    geom_hline(yintercept = 0, color = "black") +
    scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Cancer type",
      y = "Number of perturbed links\n(loss shown as negative)",
      fill = "Direction",
      title = "Gain and loss links by cancer"
    )
  
  pdf(out_pdf, width = 8, height = 5)
  print(p)
  dev.off()
}

build_pairwise_matrix <- function(df_wide, direction_label, cancer_levels) {
  mat <- matrix(
    0L,
    nrow = length(cancer_levels),
    ncol = length(cancer_levels),
    dimnames = list(cancer_levels, cancer_levels)
  )
  
  for (i in seq_along(cancer_levels)) {
    for (j in seq_along(cancer_levels)) {
      if (i == j) {
        next
      }
      ca1 <- cancer_levels[i]
      ca2 <- cancer_levels[j]
      mat[ca1, ca2] <- sum(
        df_wide[[ca1]] == direction_label &
          df_wide[[ca2]] == direction_label,
        na.rm = TRUE
      )
    }
  }
  
  mat
}

plot_pairwise_heatmap <- function(atlas_df, cancer_levels, out_pdf) {
  direction_wide <- atlas_df %>%
    dplyr::select(link, cancer, direction) %>%
    dplyr::filter(cancer %in% cancer_levels) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = cancer,
      values_from = direction,
      values_fill = NA_character_
    )
  
  gain_mat <- build_pairwise_matrix(direction_wide, "gain", cancer_levels)
  loss_mat <- build_pairwise_matrix(direction_wide, "loss", cancer_levels)
  
  plot_df <- expand.grid(
    cancer_x = cancer_levels,
    cancer_y = cancer_levels,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      x_idx = match(cancer_x, cancer_levels),
      y_idx = match(cancer_y, cancer_levels),
      gain_shared = mapply(function(y, x) gain_mat[y, x], cancer_y, cancer_x),
      loss_shared = mapply(function(y, x) loss_mat[y, x], cancer_y, cancer_x),
      fill_value = dplyr::case_when(
        x_idx > y_idx ~ gain_shared,
        x_idx < y_idx ~ -loss_shared,
        TRUE ~ NA_real_
      ),
      label = dplyr::case_when(
        x_idx > y_idx ~ as.character(gain_shared),
        x_idx < y_idx ~ as.character(loss_shared),
        TRUE ~ NA_character_
      ),
      triangle = dplyr::case_when(
        x_idx > y_idx ~ "gain",
        x_idx < y_idx ~ "loss",
        TRUE ~ "diag"
      ),
      cancer_x = factor(cancer_x, levels = cancer_levels),
      cancer_y = factor(cancer_y, levels = rev(cancer_levels))
    )
  
  max_fill <- max(c(plot_df$gain_shared, plot_df$loss_shared), na.rm = TRUE)
  
  p <- ggplot() +
    geom_tile(
      data = plot_df %>% dplyr::filter(triangle != "diag"),
      aes(x = cancer_x, y = cancer_y, fill = fill_value),
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      data = plot_df %>% dplyr::filter(triangle != "diag"),
      aes(x = cancer_x, y = cancer_y, label = label),
      size = 3.1
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-max_fill, max_fill),
      breaks = c(-max_fill, 0, max_fill),
      labels = c("Loss shared", "0", "Gain shared"),
      name = "Shared links"
    ) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank()
    ) +
    labs(
      title = "Pairwise shared gain/loss links",
      subtitle = "Lower triangle: shared gain; upper triangle: shared loss"
    )
  
  pdf(out_pdf, width = 9, height = 8)
  print(p)
  dev.off()
}

plot_gain_loss_by_cancer(atlas_all, cancers, out_gain_pdf)
plot_pairwise_heatmap(atlas_all, cancers, out_pair_pdf)

message("Saved: ", out_gain_pdf)
message("Saved: ", out_pair_pdf)
