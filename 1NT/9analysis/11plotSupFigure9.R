rm(list = ls())

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

base_family <- "Helvetica"

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
case_dir <- file.path(base_dir, "case_survival", "MMP9_SPP1")
out_file <- file.path(base_dir, "supfigure9_panel.pdf")

theme_main <- theme_bw(base_size = 11, base_family = base_family) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
        axis.title = element_text(face = "bold", family = base_family),
        axis.text = element_text(color = "#222222", family = base_family),
        plot.title = element_text(face = "bold", size = 12, family = base_family),
        plot.subtitle = element_text(size = 9.5, family = base_family),
        legend.title = element_text(face = "bold", family = base_family),
        legend.text = element_text(family = base_family),
        panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6)
    )

recurrent_tbl <- read.csv(
    file.path(base_dir, "universal_recurrent_links.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

module_annotations <- read.csv(
    file.path(base_dir, "module_results", "module_annotation_table.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

module_gene_counts <- read.csv(
    file.path(base_dir, "module_gene_counts.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

dataset_tbl <- read.csv(
    file.path(case_dir, "edge_cox_per_dataset.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Panel A: recurrence distribution
plot_a_df <- recurrent_tbl %>%
    filter(recurrent_class %in% c("recurrently_gained", "recurrently_lost")) %>%
    mutate(
        recurrence_bin = factor(recurrent_count, levels = sort(unique(recurrent_count))),
        recurrent_class = factor(
            recurrent_class,
            levels = c("recurrently_gained", "recurrently_lost"),
            labels = c("Gained", "Lost")
        )
    ) %>%
    count(recurrence_bin, recurrent_class, name = "n_links")

p_a <- ggplot(plot_a_df, aes(x = recurrence_bin, y = n_links, fill = recurrent_class)) +
    geom_col(width = 0.75) +
    geom_text(
        aes(label = n_links),
        position = position_stack(vjust = 0.5),
        size = 3.0,
        family = base_family
    ) +
    scale_fill_manual(values = c(Gained = "#b2182b", Lost = "#2166ac")) +
    labs(
        title = "Distribution of recurrent interaction rewiring",
        x = "Number of cancer types",
        y = "Number of recurrent links",
        fill = "Class"
    ) +
    theme_main +
    theme(
        legend.position = "top"
    )

# Panel B: module gene counts
plot_b_df <- module_gene_counts %>%
    group_by(direction) %>%
    slice_max(order_by = n_genes, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(direction, desc(n_genes), module_id) %>%
    mutate(
        module_id = factor(module_id, levels = rev(module_id))
    )

p_b <- ggplot(plot_b_df, aes(x = n_genes, y = module_id, fill = direction)) +
    geom_col(width = 0.72) +
    geom_text(
        aes(
            x = n_genes + max(n_genes, na.rm = TRUE) * 0.02,
            label = n_genes
        ),
        hjust = 0,
        size = 2.8,
        family = base_family
    ) +
    scale_fill_manual(values = c(gain = "#b2182b", loss = "#2166ac")) +
    labs(
        title = "Gene counts across recurrent rewiring modules",
        x = "Number of genes",
        y = "Module ID",
        fill = "Direction"
    ) +
    theme_main +
    theme(
        legend.position = "top"
    ) +
    expand_limits(x = max(plot_b_df$n_genes, na.rm = TRUE) * 1.14)

# Panel C: per-dataset survival
plot_c_df <- dataset_tbl %>%
    mutate(
        dataset_label = paste(cancer, dataset, sep = " | "),
        dataset_label = factor(dataset_label, levels = dataset_label[order(hr)]),
        direction = if_else(hr >= 1, "risk", "protective")
    )

p_c <- ggplot(plot_c_df, aes(x = hr, y = dataset_label, color = direction)) +
    geom_point(size = 1.9) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.15, linewidth = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#8c8c8c", linewidth = 0.5) +
    scale_x_log10() +
    scale_color_manual(values = c(risk = "#b2182b", protective = "#2166ac")) +
    labs(
        title = "MMP9-SPP1 survival effect across individual datasets",
        x = "Hazard ratio",
        y = "Cancer | Dataset",
        color = "Direction"
    ) +
    theme_main +
    theme(
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 4.0, family = base_family),
        legend.position = "top"
    )

p_left <- p_a / p_b

p_final <- (p_left | p_c) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(face = "bold", size = 14, family = base_family)
    )

ggsave(out_file, p_final, width = 15, height = 10)

cat("Saved ->", out_file, "\n")
