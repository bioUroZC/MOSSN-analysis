rm(list = ls())

library(ggplot2)
library(dplyr)
library(data.table)
library(cluster)
library(tidyr)
library(patchwork)

base_dir <- "/proj/c.zihao/work1/1NT/8dataDriven/1distance"
out_dir  <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_df <- read.csv(file.path(base_dir, "results_Cluster.csv"), stringsAsFactors = FALSE)

cancer_levels   <- c("LUAD")
fraction_levels <- c("5%", "10%", "15%", "20%")

method_levels <- c("MOSSN", "SSN", "SWEET", "LIONESS")

cluster_df <- cluster_df %>%
    filter(method %in% method_levels, cancer %in% cancer_levels) %>%
    mutate(
        cancer           = factor(cancer, levels = cancer_levels),
        method           = factor(method, levels = method_levels),
        feature_fraction = factor(feature_fraction, levels = fraction_levels)
    )

heatmap_theme <- theme_minimal(base_size = 12) +
    theme(
        panel.grid      = element_blank(),
        axis.title      = element_text(size = 12, face = "bold"),
        axis.text.x     = element_text(angle = 30, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y     = element_text(face = "bold"),
        legend.title    = element_text(face = "bold"),
        legend.key.height = unit(42, "pt")
    )

auc_summary <- cluster_df %>%
    group_by(feature_fraction, method) %>%
    summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop")

accuracy_summary <- cluster_df %>%
    group_by(feature_fraction, method) %>%
    summarise(mean_accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

auc_limits <- range(auc_summary$mean_auc, na.rm = TRUE)
acc_limits <- range(accuracy_summary$mean_accuracy, na.rm = TRUE)

p_auc_method_fraction <- ggplot(auc_summary, aes(x = method, y = feature_fraction, fill = mean_auc)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", mean_auc)), size = 3.0) +
    scale_fill_gradientn(
        colours = c("#fff8e8", "#fdb366", "#d94701"),
        limits = auc_limits,
        breaks = pretty(auc_limits, n = 4)
    ) +
    labs(x = "Method", y = "Feature fraction", fill = "Mean AUC") +
    coord_fixed() +
    heatmap_theme

ggsave(file.path(out_dir, "cluster_auc_method_fraction_heatmap.pdf"),
       p_auc_method_fraction, width = 7.0, height = 4.0)

p_accuracy_method_fraction <- ggplot(accuracy_summary, aes(x = method, y = feature_fraction, fill = mean_accuracy)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", mean_accuracy)), size = 3.0) +
    scale_fill_gradientn(
        colours = c("#f4fbf2", "#8fd19e", "#0b6e4f"),
        limits = acc_limits,
        breaks = pretty(acc_limits, n = 4)
    ) +
    labs(x = "Method", y = "Feature fraction", fill = "Mean accuracy") +
    coord_fixed() +
    heatmap_theme

ggsave(file.path(out_dir, "cluster_accuracy_method_fraction_heatmap.pdf"),
       p_accuracy_method_fraction, width = 7.0, height = 4.0)

write.csv(auc_summary,      file.path(out_dir, "cluster_auc_overall_summary.csv"),       row.names = FALSE)
write.csv(accuracy_summary, file.path(out_dir, "cluster_accuracy_overall_summary.csv"),  row.names = FALSE)

# ── PCA plots ────────────────────────────────────────────────────────────────
metadata_path <- "/proj/c.zihao/work1/1NT/1data/TCGA/metadata.csv"
mat_dir       <- "/proj/c.zihao/work1/1NT/8dataDriven"
pca_cancer    <- "LUAD"
pca_top_frac  <- 0.20
pca_seed      <- 123

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
cancer_samples <- metadata$Sample[metadata$Type == pca_cancer]

sep_metrics_list <- list()
pca_plot_list    <- list()

pca_theme <- theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "bottom"
    )

for (method in method_levels) {
    mat_path <- file.path(mat_dir, method, "merged_matrix.csv")
    if (!file.exists(mat_path)) {
        cat("Skipping", method, ": merged_matrix.csv not found\n")
        next
    }

    mat    <- fread(mat_path) |> as.data.frame()
    id_col <- intersect(c("Interaction", "V1", "Unnamed: 0"), colnames(mat))[1]
    rownames(mat) <- mat[[id_col]]
    mat[[id_col]] <- NULL
    mat <- abs(mat)

    cols <- intersect(cancer_samples, colnames(mat))
    if (length(cols) < 4) next
    sub_mat <- mat[, cols, drop = FALSE]
    sub_mat[is.na(sub_mat)] <- 0

    cv   <- apply(sub_mat, 1, sd) / abs(rowMeans(sub_mat))
    cv   <- cv[is.finite(cv)]
    top_n <- floor(length(cv) * pca_top_frac)
    keep  <- names(sort(cv, decreasing = TRUE))[seq_len(top_n)]

    x_raw <- t(as.matrix(sub_mat[keep, , drop = FALSE]))
    x_raw <- x_raw[, apply(x_raw, 2, sd) > 0, drop = FALSE]
    x_sc  <- scale(x_raw)

    set.seed(pca_seed)
    pca_res <- prcomp(x_sc, center = FALSE, scale. = FALSE)
    var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

    df_pca <- data.frame(
        pc1  = pca_res$x[, 1],
        pc2  = pca_res$x[, 2],
        type = ifelse(grepl("11A$", rownames(x_sc)), "Normal", "Tumor")
    )

    p_pca <- ggplot(df_pca, aes(pc1, pc2, color = type)) +
        geom_point(size = 2, alpha = 0.8) +
        stat_ellipse(
            type = "norm", linewidth = 0.6,
            linetype = "dashed", level = 0.90
        ) +
        scale_color_manual(
            values = c(Tumor = "#d62728", Normal = "#1f77b4")
        ) +
        labs(
            title = sprintf("%s - %s", pca_cancer, method),
            x     = sprintf("PC1 (%.1f%%)", var_exp[1]),
            y     = sprintf("PC2 (%.1f%%)", var_exp[2]),
            color = NULL
        ) +
        pca_theme

    fname <- sprintf("pca_%s_%s.pdf", pca_cancer, method)
    ggsave(file.path(out_dir, fname), p_pca, width = 6, height = 6)
    pca_plot_list[[method]] <- p_pca
    cat("PCA saved:", method, "\n")

    # Separation metrics in PC1-PC2 space
    pc12   <- as.matrix(df_pca[, c("pc1", "pc2")])
    labels <- ifelse(df_pca$type == "Tumor", 1L, 0L)

    sil     <- silhouette(labels, dist(pc12))
    sil_all <- mean(sil[, "sil_width"])
    sil_t   <- mean(sil[labels == 1, "sil_width"])
    sil_n   <- mean(sil[labels == 0, "sil_width"])

    var_exp_12 <- sum(var_exp[1:2])

    sep_metrics_list[[method]] <- data.frame(
        method     = method,
        sil_tumor  = round(sil_t,       4),
        sil_normal = round(sil_n,       4),
        sil_mean   = round(sil_all,     4),
        var_exp_12 = round(var_exp_12,  2)
    )
    cat(sprintf(
        "  sil_tumor=%.3f  sil_normal=%.3f  var_exp12=%.1f%%\n",
        sil_t, sil_n, var_exp_12
    ))
}

sep_metrics <- do.call(rbind, sep_metrics_list)
sep_metrics$method <- factor(sep_metrics$method, levels = method_levels)
write.csv(sep_metrics,
          file.path(out_dir, "pca_separation_metrics.csv"),
          row.names = FALSE)

# Silhouette by group (Normal vs Tumor)
sil_long <- tidyr::pivot_longer(
    sep_metrics,
    cols      = c("sil_tumor", "sil_normal"),
    names_to  = "group",
    values_to = "silhouette"
)
sil_long$group <- factor(
    sil_long$group,
    levels = c("sil_normal", "sil_tumor"),
    labels = c("Normal", "Tumor")
)

bar_colors <- c(Normal = "#1f77b4", Tumor = "#d62728")

p_sil <- ggplot(sil_long,
                aes(x = method, y = silhouette, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    scale_fill_manual(values = bar_colors) +
    labs(
        x     = NULL, y = "Silhouette score", fill = NULL,
        title = sprintf(
            "Group silhouette in PC1-PC2 space (%s)", pca_cancer
        )
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5),
        legend.position    = "bottom",
        axis.text.x        = element_text(angle = 30, hjust = 1)
    )

ggsave(file.path(out_dir, "pca_silhouette_by_group.pdf"),
       p_sil, width = 5, height = 4.5)

# Variance explained by PC1+PC2
p_var <- ggplot(sep_metrics,
                aes(x = method, y = var_exp_12, fill = method)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(
        MOSSN = "#2166ac", SSN = "#f4a582",
        SWEET = "#92c5de", LIONESS = "#d6604d"
    )) +
    labs(
        x = NULL, y = "Variance explained (%)",
        title = sprintf("PC1+PC2 variance explained (%s)", pca_cancer)
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5),
        legend.position    = "none",
        axis.text.x        = element_text(angle = 30, hjust = 1)
    )

ggsave(file.path(out_dir, "pca_variance_explained.pdf"),
       p_var, width = 4.5, height = 4.5)

# ── Combined PCA figure: 4 methods in a 2x2 layout ───────────────────────────
pca_panel_titles <- c(
    MOSSN   = "MOSSN-driven",
    SSN     = "SSN",
    SWEET   = "SWEET",
    LIONESS = "LIONESS"
)

pca_panels_2x2 <- lapply(method_levels, function(method) {
    pca_plot_list[[method]] +
        labs(title = pca_panel_titles[[method]]) +
        theme(
            plot.title      = element_text(face = "bold", hjust = 0.5),
            legend.position = "bottom"
        )
})
names(pca_panels_2x2) <- method_levels

p_pca_2x2 <- wrap_plots(pca_panels_2x2, ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(face = "bold", size = 14),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "pca_LUAD_4methods_2x2.pdf"),
       p_pca_2x2, width = 12, height = 10)

# ── Figure 1: MOSSN + SSN + metrics (4 panels, one row) ──────────────────────
# Shorter titles for the combined figure
p_var_fig <- p_var +
    labs(title = "PC1+PC2 variance explained")
p_sil_fig <- p_sil +
    labs(title = "Silhouette by group")

p_fig1 <- (
    pca_plot_list[["MOSSN"]] +
    pca_plot_list[["SSN"]] +
    p_var_fig +
    p_sil_fig
) +
    plot_layout(nrow = 1, widths = c(1.3, 1.3, 1, 1.3)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, "fig1_main.pdf"),
       p_fig1, width = 20, height = 6)

# ── Figure 2: SWEET + LIONESS (2 panels, one row) ────────────────────────────
p_fig2 <- (
    pca_plot_list[["SWEET"]] +
    pca_plot_list[["LIONESS"]]
) +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, "fig2_supp.pdf"),
       p_fig2, width = 10, height = 5.5)
