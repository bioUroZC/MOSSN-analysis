rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
    library(ggplot2)
})

base_family <- "Helvetica"

RESULTS_DIR <- paste0(PROJ_ROOT, "/4coupled/results")
IN_DIR <- file.path(RESULTS_DIR, "survival_stratification")
OUT_DIR <- file.path(IN_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_CANCERS <- c("BLCA", "LIHC", "LUAD", "SARC", "STAD")
METHODS <- c(
    "MOSSN_EXP", "MOSSN_MET", "MOSSN_CNV",
    "MOSSN_NoCross", "MOSSN_Restart",
    "MOSSN_DirectNoDyn", "MOSSN_MultiLayer", "MOSSN_Direct"
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
    "MOSSN_EXP"         = "#4daf4a",
    "MOSSN_MET"         = "#a65628",
    "MOSSN_CNV"         = "#f781bf",
    "MOSSN_NoCross"     = "#999999",
    "MOSSN_Restart"     = "#377eb8",
    "MOSSN_DirectNoDyn" = "#ff7f00",
    "MOSSN_MultiLayer"  = "#984ea3",
    "MOSSN_Direct"      = "#e41a1c"
)

method_cols_labeled <- setNames(unname(method_cols[METHODS]), method_labels[METHODS])

all_df <- read.csv(file.path(IN_DIR, "balanced_results.csv"), stringsAsFactors = FALSE)
sub_df <- all_df[
    all_df$Cancer %in% TARGET_CANCERS &
    all_df$Method %in% METHODS &
    all_df$FeatureFrac %in% FEATURE_FRACS,
]

sub_df$FeaturePct <- as.numeric(sub("%", "", sub_df$FeatureFrac))

avg_df <- aggregate(
    NegLog10_LogRank ~ Method + FeaturePct,
    data = sub_df,
    FUN = mean,
    na.rm = TRUE
)
avg_df$Method <- factor(avg_df$Method, levels = METHODS)
avg_df$MethodLabel <- factor(method_labels[avg_df$Method], levels = method_labels[METHODS])
avg_df$Highlight <- ifelse(avg_df$Method == "MOSSN_Direct", "highlight", "background")

plot_theme <- theme_bw(base_size = 12, base_family = base_family) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
        axis.title.x = element_text(face = "bold", family = base_family),
        axis.title.y = element_text(face = "bold", family = base_family),
        axis.text = element_text(color = "#222222", family = base_family),
        plot.title = element_text(face = "bold", size = 14, family = base_family),
        plot.subtitle = element_text(size = 10, color = "#444444", family = base_family),
        legend.position = "top",
        legend.text = element_text(size = 9, family = base_family),
        legend.spacing.x = unit(4, "pt"),
        legend.key.width = unit(12, "pt"),
        legend.key.height = unit(12, "pt"),
        panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6),
        plot.margin = margin(8, 10, 8, 8)
    )

p <- ggplot() +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        linewidth = 0.5,
        color = "#8c8c8c"
    ) +
    geom_line(
        data = avg_df[avg_df$Highlight == "background", ],
        aes(x = FeaturePct, y = NegLog10_LogRank, color = MethodLabel, group = MethodLabel),
        linewidth = 0.85,
        alpha = 0.95
    ) +
    geom_point(
        data = avg_df[avg_df$Highlight == "background", ],
        aes(x = FeaturePct, y = NegLog10_LogRank, color = MethodLabel),
        size = 2.4,
        alpha = 0.98
    ) +
    geom_line(
        data = avg_df[avg_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = NegLog10_LogRank, color = MethodLabel, group = MethodLabel),
        linewidth = 1.6
    ) +
    geom_point(
        data = avg_df[avg_df$Highlight == "highlight", ],
        aes(x = FeaturePct, y = NegLog10_LogRank, color = MethodLabel),
        size = 3.2
    ) +
    scale_x_continuous(
        breaks = c(5, 10, 15, 20),
        expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = method_cols_labeled, drop = FALSE) +
    labs(
        title = "Log-rank significance across feature fractions",
        subtitle = "Mean -log10(P) across BLCA, LIHC, LUAD, SARC, and STAD",
        x = "Feature percentage",
        y = "-log10(P)",
        color = NULL
    ) +
    plot_theme +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.2)))

out_file <- file.path(OUT_DIR, "balanced_mean_logp_lines_selected_cancers.pdf")
ggsave(out_file, p, width = 6, height = 6)
cat("Saved plot ->", out_file, "\n")
