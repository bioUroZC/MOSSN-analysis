rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

base_family <- "Helvetica"

RESULTS_DIR <- paste0(PROJ_ROOT, "/4coupled/results")
IN_DIR <- file.path(RESULTS_DIR, "survival_stratification")
OUT_DIR <- file.path(IN_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_CANCERS <- c("CRC")
METHODS <- c(
    "MOSSN_EXP", "MOSSN_MET", "MOSSN_CNV",
    "MOSSN_NoCross", "MOSSN_Restart",
    "MOSSN_DirectNoDyn", "MOSSN_MultiLayer", "MOSSN_Direct"
)
legacy_to_new <- c(
    "EXP_single" = "MOSSN_EXP",
    "MET_single" = "MOSSN_MET",
    "CNV_single" = "MOSSN_CNV",
    "MUL_noCross" = "MOSSN_NoCross",
    "MUL_full" = "MOSSN_Restart",
    "MUL_direct" = "MOSSN_Direct",
    "MUL_direct_fixed" = "MOSSN_DirectNoDyn",
    "MUL_multilayer" = "MOSSN_MultiLayer"
)
FEATURE_FRACS <- c("5%", "10%", "15%", "20%")

method_labels <- c(
    "MOSSN_EXP" = "MOSSN-EXP",
    "MOSSN_MET" = "MOSSN-MET",
    "MOSSN_CNV" = "MOSSN-CNV",
    "MOSSN_NoCross" = "MOSSN-NoCross",
    "MOSSN_Restart" = "MOSSN-Restart",
    "MOSSN_DirectNoDyn" = "MOSSN-DirectNoDyn",
    "MOSSN_MultiLayer" = "MOSSN-MultiLayer",
    "MOSSN_Direct" = "MOSSN-Direct"
)

method_cols <- c(
    "MOSSN_EXP"         = "#4daf4a",   # green
    "MOSSN_MET"         = "#a65628",   # brown
    "MOSSN_CNV"         = "#f781bf",   # pink
    "MOSSN_NoCross"     = "#999999",   # gray
    "MOSSN_Restart"     = "#377eb8",   # blue
    "MOSSN_DirectNoDyn" = "#ff7f00",   # orange
    "MOSSN_MultiLayer"  = "#984ea3",   # purple
    "MOSSN_Direct"      = "#e41a1c"    # red — main method
)

method_cols_labeled <- setNames(unname(method_cols[METHODS]), method_labels[METHODS])

all_df <- read.csv(
    file.path(IN_DIR, "balanced_results.csv"),
    stringsAsFactors = FALSE
)
all_df$Method <- ifelse(all_df$Method %in% names(legacy_to_new),
                        legacy_to_new[all_df$Method],
                        all_df$Method)

plot_theme <- theme_bw(base_size = 12, base_family = base_family) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
        axis.title.x = element_text(face = "bold", family = base_family),
        axis.title.y = element_text(face = "bold", family = base_family),
        axis.text = element_text(color = "#222222", family = base_family),
        plot.title = element_text(face = "bold", size = 14, family = base_family),
        legend.position = "top",
        legend.text = element_text(size = 9, family = base_family),
        legend.spacing.x = unit(4, "pt"),
        legend.key.width = unit(14, "pt"),
        panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6),
        plot.margin = margin(8, 10, 8, 8)
    )

plot_guides <- guides(
    color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.2)
    )
)

# Filter to target cancers and compute mean across cancers per Method x FeatureFrac
sub_df <- all_df[
    all_df$Cancer %in% TARGET_CANCERS &
    all_df$Method %in% METHODS &
    all_df$FeatureFrac %in% FEATURE_FRACS,
]
sub_df$FeaturePct <- as.numeric(sub("%", "", sub_df$FeatureFrac))
sub_df$NegLog10P  <- -log10(sub_df$LogRank_p)

avg_df <- aggregate(
    cbind(NegLog10P, Cox_HR_Poor_vs_Good) ~ Method + FeaturePct,
    data = sub_df,
    FUN  = mean,
    na.rm = TRUE
)
avg_df$Method      <- factor(avg_df$Method, levels = METHODS)
avg_df$MethodLabel <- factor(method_labels[avg_df$Method], levels = method_labels[METHODS])
avg_df$Highlight   <- ifelse(avg_df$Method == "MOSSN_Direct", "highlight", "background")

# Log-rank panel
logrank_df <- avg_df[, c("Method", "MethodLabel", "FeaturePct", "Highlight", "NegLog10P")]
colnames(logrank_df)[5] <- "Value"

p_logrank <- ggplot() +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        linewidth = 0.5,
        color = "#8c8c8c"
    ) +
    geom_line(
        data = logrank_df[logrank_df$Highlight == "background", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel, group = MethodLabel),
        linewidth = 0.8,
        alpha = 0.9
    ) +
    geom_point(
        data = logrank_df[logrank_df$Highlight == "background", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel),
        size = 2.3,
        alpha = 0.95
    ) +
    geom_line(
        data = logrank_df[logrank_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel, group = MethodLabel),
        linewidth = 1.5
    ) +
    geom_point(
        data = logrank_df[logrank_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel),
        size = 3.2
    ) +
    scale_x_continuous(
        breaks = sort(unique(logrank_df$FeaturePct)),
        expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = method_cols_labeled, drop = FALSE) +
    labs(
        title = "Log-rank significance (mean across cancers)",
        x = "Feature percentage",
        y = "-log10(P)",
        color = NULL
    ) +
    plot_theme +
    plot_guides

# Hazard ratio panel
hr_df <- avg_df[, c("Method", "MethodLabel", "FeaturePct", "Highlight", "Cox_HR_Poor_vs_Good")]
colnames(hr_df)[5] <- "Value"

p_hr <- ggplot() +
    geom_hline(
        yintercept = 1,
        linetype = "dashed",
        linewidth = 0.5,
        color = "#8c8c8c"
    ) +
    geom_line(
        data = hr_df[hr_df$Highlight == "background", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel, group = MethodLabel),
        linewidth = 0.8,
        alpha = 0.9
    ) +
    geom_point(
        data = hr_df[hr_df$Highlight == "background", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel),
        size = 2.3,
        alpha = 0.95
    ) +
    geom_line(
        data = hr_df[hr_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel, group = MethodLabel),
        linewidth = 1.5
    ) +
    geom_point(
        data = hr_df[hr_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = Value, color = MethodLabel),
        size = 3.2
    ) +
    scale_x_continuous(
        breaks = sort(unique(hr_df$FeaturePct)),
        expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = method_cols_labeled, drop = FALSE) +
    labs(
        title = "Hazard ratio (mean across cancers)",
        x = "Feature percentage",
        y = "Hazard ratio",
        color = NULL
    ) +
    plot_theme +
    plot_guides

p <- p_logrank + p_hr + plot_layout(ncol = 2)

out_file <- file.path(OUT_DIR, "balanced_mean_survival_two_panel.pdf")
ggsave(out_file, p, width = 12, height = 5.8)
cat("Saved plot ->", out_file, "\n")
