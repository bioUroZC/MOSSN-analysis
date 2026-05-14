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
out_file <- file.path(base_dir, "figure5_panel.pdf")

cancers <- c(
    "BLCA", "BRCA", "CRC", "ESCA", "HNSC",
    "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD"
)

theme_main <- theme_bw(base_size = 11, base_family = base_family) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
        axis.title = element_text(face = "bold", family = base_family),
        axis.text = element_text(color = "#222222", family = base_family),
        plot.title = element_text(face = "bold", size = 12, family = base_family),
        plot.subtitle = element_text(size = 9.5, family = base_family),
        legend.position = "right",
        legend.title = element_text(face = "bold", family = base_family),
        legend.text = element_text(family = base_family),
        panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6)
    )

atlas_all <- read.csv(
    file.path(base_dir, "atlas_all.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if ("" %in% names(atlas_all)) {
    atlas_all <- atlas_all[, names(atlas_all) != "", drop = FALSE]
}

module_edges <- read.csv(
    file.path(case_dir, "edge_case_module_edges.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

module_annotations <- read.csv(
    file.path(base_dir, "module_results", "module_annotation_table.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

case_summary <- read.csv(
    file.path(case_dir, "edge_case_summary.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

meta_summary <- read.csv(
    file.path(case_dir, "edge_cox_meta_summary.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

priority_tbl <- read.csv(
    file.path(base_dir, "edge_prioritization_candidates.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Panel A
plot_a_df <- atlas_all %>%
    filter(direction %in% c("gain", "loss"), cancer %in% cancers) %>%
    count(cancer, direction, name = "n_links") %>%
    complete(
        cancer = cancers,
        direction = c("gain", "loss"),
        fill = list(n_links = 0L)
    ) %>%
    mutate(
        cancer = factor(cancer, levels = cancers),
        direction = factor(direction, levels = c("gain", "loss")),
        y = if_else(direction == "loss", -n_links, n_links)
    )

offset_a <- max(plot_a_df$n_links, na.rm = TRUE) * 0.02

p_a <- ggplot(plot_a_df, aes(x = cancer, y = y, fill = direction)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_text(
        aes(
            label = n_links,
            y = ifelse(direction == "gain", y + offset_a, y - offset_a)
        ),
        size = 2.8,
        family = base_family
    ) +
    scale_fill_manual(values = c(gain = "#b2182b", loss = "#2166ac")) +
    labs(
        title = "Pan-cancer gain/loss burden",
        x = "Cancer type",
        y = "Number of perturbed links\n(loss shown as negative)"
    ) +
    theme_main +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

# Panel B
direction_wide <- atlas_all %>%
    select(link, cancer, direction) %>%
    filter(cancer %in% cancers) %>%
    distinct() %>%
    pivot_wider(
        names_from = cancer,
        values_from = direction,
        values_fill = NA_character_
    )

gain_mat <- matrix(0L, nrow = length(cancers), ncol = length(cancers), dimnames = list(cancers, cancers))
loss_mat <- matrix(0L, nrow = length(cancers), ncol = length(cancers), dimnames = list(cancers, cancers))

for (i in seq_along(cancers)) {
    for (j in seq_along(cancers)) {
        if (i == j) {
            next
        }
        ca1 <- cancers[i]
        ca2 <- cancers[j]
        gain_mat[ca1, ca2] <- sum(
            direction_wide[[ca1]] == "gain" & direction_wide[[ca2]] == "gain",
            na.rm = TRUE
        )
        loss_mat[ca1, ca2] <- sum(
            direction_wide[[ca1]] == "loss" & direction_wide[[ca2]] == "loss",
            na.rm = TRUE
        )
    }
}

plot_b_df <- expand.grid(
    cancer_x = cancers,
    cancer_y = cancers,
    stringsAsFactors = FALSE
) %>%
    mutate(
        x_idx = match(cancer_x, cancers),
        y_idx = match(cancer_y, cancers),
        gain_shared = mapply(function(y, x) gain_mat[y, x], cancer_y, cancer_x),
        loss_shared = mapply(function(y, x) loss_mat[y, x], cancer_y, cancer_x),
        fill_value = case_when(
            x_idx > y_idx ~ gain_shared,
            x_idx < y_idx ~ -loss_shared,
            TRUE ~ NA_real_
        ),
        label = case_when(
            x_idx > y_idx ~ as.character(gain_shared),
            x_idx < y_idx ~ as.character(loss_shared),
            TRUE ~ NA_character_
        ),
        triangle = case_when(
            x_idx > y_idx ~ "gain",
            x_idx < y_idx ~ "loss",
            TRUE ~ "diag"
        ),
        cancer_x = factor(cancer_x, levels = cancers),
        cancer_y = factor(cancer_y, levels = rev(cancers))
    )

max_fill_b <- max(c(plot_b_df$gain_shared, plot_b_df$loss_shared), na.rm = TRUE)

p_b <- ggplot() +
    geom_tile(
        data = plot_b_df %>% filter(triangle != "diag"),
        aes(x = cancer_x, y = cancer_y, fill = fill_value),
        color = "white",
        linewidth = 0.35
    ) +
    geom_text(
        data = plot_b_df %>% filter(triangle != "diag"),
        aes(x = cancer_x, y = cancer_y, label = label),
        size = 2.6,
        family = base_family
    ) +
    scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-max_fill_b, max_fill_b),
        breaks = c(-max_fill_b, 0, max_fill_b),
        labels = c("Shared loss", "0", "Shared gain"),
        name = "Shared links"
    ) +
    coord_fixed() +
    labs(
        title = "Shared recurrent rewiring across cancers",
        x = NULL,
        y = NULL
    ) +
    theme_main +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()
    )

# Panel C
plot_c_df <- module_annotations %>%
    filter(annotation_status == "ok") %>%
    group_by(direction) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(direction, module_rank) %>%
    mutate(
        module_id = factor(module_id, levels = rev(module_id))
    )

p_c <- ggplot(plot_c_df, aes(x = n_edges, y = module_id, fill = direction)) +
    geom_col(width = 0.72) +
    geom_text(
        aes(
            x = pmax(n_edges * 0.03, 0.8),
            label = concise_label
        ),
        hjust = 0,
        size = 2.8,
        family = base_family,
        color = "black"
    ) +
    scale_fill_manual(values = c(gain = "#b2182b", loss = "#2166ac")) +
    labs(
        title = "Top annotated rewiring modules",
        x = "Number of recurrent edges",
        y = "Module ID"
    ) +
    theme_main +
    theme(
        legend.position = "top"
    )

# Panel D
label_tbl_d <- bind_rows(
    priority_tbl %>% filter(is_target),
    priority_tbl %>% filter(is_annotated) %>% slice_head(n = 6)
) %>%
    distinct(link, .keep_all = TRUE)

priority_tbl$highlight_group <- case_when(
    priority_tbl$is_target ~ "Target edge",
    priority_tbl$link %in% label_tbl_d$link ~ "Prioritized edges",
    TRUE ~ "Background"
)

p_d <- ggplot(priority_tbl, aes(x = recurrent_count, y = weight)) +
    geom_point(
        data = priority_tbl %>% filter(highlight_group == "Background"),
        color = "grey82",
        alpha = 0.7,
        size = 1.8
    ) +
    geom_point(
        data = priority_tbl %>% filter(highlight_group == "Prioritized edges"),
        aes(color = direction),
        alpha = 0.95,
        size = 3.0
    ) +
    geom_point(
        data = priority_tbl %>% filter(is_target),
        color = "black",
        fill = "#ffd54f",
        shape = 21,
        size = 4.0,
        stroke = 1.0
    ) +
    geom_text(
        data = label_tbl_d,
        aes(label = link),
        size = 2.8,
        nudge_y = 0.18,
        check_overlap = TRUE,
        color = "black",
        family = base_family,
        show.legend = FALSE
    ) +
    scale_color_manual(values = c(gain = "#b2182b", loss = "#2166ac")) +
    labs(
        title = "Prioritization of recurrent interactions",
        x = "Recurrent cancer count",
        y = "Recurrence × |Δ weight| (tumor vs normal)"
    ) +
    theme_main +
    theme(
        legend.position = "none"
    )

# Panel E
plot_e_df <- module_edges %>%
    slice_head(n = 12) %>%
    mutate(
        edge_label = factor(link, levels = rev(link)),
        fill_group = if_else(is_selected_edge, "selected", "context")
    )

p_e <- ggplot(plot_e_df, aes(x = weight, y = edge_label, fill = fill_group)) +
    geom_col(width = 0.72) +
    scale_fill_manual(
        values = c(selected = "#b2182b", context = "#bdbdbd"),
        labels = c(selected = "MMP9-SPP1", context = "Other module edges")
    ) +
    labs(
        title = "MMP9-SPP1 in recurrent gain module",
        subtitle = paste(case_summary$module_id[1], "|", case_summary$concise_label[1]),
        x = "Module edge weight",
        y = "Edges"
    ) +
    theme_main +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )

# Panel F
plot_f_df <- meta_summary %>%
    filter(scope != "all_datasets") %>%
    mutate(
        scope = factor(scope, levels = scope[order(hr_meta)]),
        direction = if_else(hr_meta >= 1, "risk", "protective")
    )

plot_f_all <- meta_summary %>%
    filter(scope == "all_datasets") %>%
    mutate(
        scope = "Pan-cancer",
        direction = if_else(hr_meta >= 1, "risk", "protective")
    )

p_f <- ggplot(plot_f_df, aes(x = hr_meta, y = scope)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#8c8c8c", linewidth = 0.5) +
    geom_errorbarh(
        aes(xmin = conf_low, xmax = conf_high, color = direction),
        height = 0.18,
        linewidth = 0.6
    ) +
    geom_point(aes(color = direction), size = 2.3) +
    geom_errorbarh(
        data = plot_f_all,
        aes(xmin = conf_low, xmax = conf_high),
        height = 0.24,
        linewidth = 0.9,
        color = "#b2182b"
    ) +
    geom_point(
        data = plot_f_all,
        size = 3.2,
        color = "#b2182b"
    ) +
    geom_text(
        data = plot_f_all,
        aes(x = conf_high, y = 0.4, label = "Pan-cancer"),
        inherit.aes = FALSE,
        hjust = -0.1,
        size = 3.0,
        family = base_family
    ) +
    scale_color_manual(values = c(risk = "#b2182b", protective = "#2166ac")) +
    labs(
        title = "Pan-cancer prognostic effect of MMP9-SPP1",
        x = "Hazard ratio",
        y = "Cancer type"
    ) +
    theme_main +
    theme(
        legend.position = "none",
        plot.margin = margin(5.5, 30, 5.5, 5.5)
    )

p_final <- (p_a + p_b + p_c) / (p_d + p_e + p_f) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(face = "bold", size = 14, family = base_family)
    )

ggsave(out_file, p_final, width = 18, height = 10.5)

cat("Saved ->", out_file, "\n")
