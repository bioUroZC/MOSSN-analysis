rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

RESULTS_DIR <- paste0(PROJ_ROOT, "/4coupled/results")
IN_FILE <- file.path(RESULTS_DIR, "survival_stratification", "balanced_results.csv")
OUT_DIR <- file.path(RESULTS_DIR, "survival_stratification", "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FEATURE_FRAC <- Sys.getenv("FEATURE_FRAC", unset = "5%")

METHODS <- c(
    "MOSSN_Direct", "MOSSN_DirectNoDyn", "MOSSN_MultiLayer",
     "MOSSN_Restart", "MOSSN_NoCross", "MOSSN_EXP", "MOSSN_MET", "MOSSN_CNV"
)

TARGET_CANCERS <- c("LUAD", "LIHC",  "STAD", "BLCA", "SARC")

method_labels <- c(
    "MOSSN_Direct" = "MOSSN-Direct",
    "MOSSN_DirectNoDyn" = "MOSSN-DirectNoDyn",
    "MOSSN_MultiLayer" = "MOSSN-MultiLayer",
    "MOSSN_Restart" = "MOSSN-Restart",
    "MOSSN_NoCross" = "MOSSN-NoCross",
    "MOSSN_EXP" = "MOSSN-EXP",
    "MOSSN_MET" = "MOSSN-MET",
    "MOSSN_CNV" = "MOSSN-CNV"
)

df <- read.csv(IN_FILE, stringsAsFactors = FALSE)
df <- df[df$Method %in% METHODS, ]
df <- df[df$FeatureFrac == FEATURE_FRAC, ]
df <- df[df$Cancer %in% TARGET_CANCERS, ]

direct_df <- df[df$Method == "MOSSN_Direct", c("Cancer", "NegLog10_LogRank")]
colnames(direct_df)[2] <- "DirectScore"

cancer_order <- direct_df$Cancer[order(direct_df$DirectScore, decreasing = TRUE)]
summary_df <- aggregate(
    NegLog10_LogRank ~ Method,
    data = df,
    FUN = mean,
    na.rm = TRUE
)
summary_df$Cancer <- "Mean"
summary_df$FeatureFrac <- FEATURE_FRAC
summary_df$CellLabel <- sprintf("%.2f", summary_df$NegLog10_LogRank)
summary_df$TextColor <- ifelse(summary_df$NegLog10_LogRank >= 2.5, "white", "#111111")

df$Cancer <- factor(df$Cancer, levels = c(cancer_order, "Mean"))
df$Method <- factor(df$Method, levels = rev(METHODS))
df$CellLabel <- sprintf("%.2f", df$NegLog10_LogRank)
df$TextColor <- ifelse(df$NegLog10_LogRank >= 2.5, "white", "#111111")

summary_df$Cancer <- factor(summary_df$Cancer, levels = levels(df$Cancer))
summary_df$Method <- factor(summary_df$Method, levels = rev(METHODS))

plot_df <- rbind(
    df[, c("Cancer", "Method", "NegLog10_LogRank", "CellLabel", "TextColor")],
    summary_df[, c("Cancer", "Method", "NegLog10_LogRank", "CellLabel", "TextColor")]
)

plot_df$Cancer <- factor(plot_df$Cancer, levels = levels(df$Cancer))
plot_df$Method <- factor(plot_df$Method, levels = rev(METHODS))

mean_col_x <- length(levels(plot_df$Cancer))
main_max <- max(df$NegLog10_LogRank, na.rm = TRUE)
fill_upper <- max(2.5, ceiling(main_max * 10) / 10)
n_cols <- length(levels(plot_df$Cancer))
n_rows <- length(levels(plot_df$Method))
panel_ratio <- n_cols / n_rows

plot_theme <- theme_minimal(base_size = 11.5) +
    theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, color = "#222222", size = 10.5),
        axis.text.y = element_text(color = "#222222", size = 10.5),
        axis.ticks = element_blank(),
        legend.position = "none"
    )

p <- ggplot(plot_df, aes(x = Cancer, y = Method, fill = NegLog10_LogRank)) +
    annotate(
        "rect",
        xmin = mean_col_x - 0.5, xmax = mean_col_x + 0.5,
        ymin = 0.5, ymax = length(METHODS) + 0.5,
        fill = "#f4f1ea", color = NA
    ) +
    geom_tile(color = "#f7f7f7", linewidth = 0.9, width = 0.92, height = 0.78) +
    geom_vline(xintercept = mean_col_x - 0.5, color = "#7f7f7f", linewidth = 0.55) +
    geom_text(
        aes(label = CellLabel, color = TextColor),
        size = 3.0,
        fontface = "plain"
    ) +
    scale_color_identity() +
    scale_y_discrete(labels = method_labels[rev(METHODS)]) +
    scale_fill_gradientn(
        colours = c("#fbf8f3", "#dceee8", "#95c7c4", "#3f8f9b", "#0f4c5c"),
        values = scales::rescale(c(0, 0.8, 1.3, 2.0, fill_upper)),
        limits = c(0, fill_upper),
        name = expression(-log[10](P))
    ) +
    labs(title = NULL) +
    coord_fixed(ratio = 0.72) +
    plot_theme

frac_tag <- sub("%", "pct", FEATURE_FRAC)
out_file <- file.path(OUT_DIR, paste0("balanced_survival_heatmap_", frac_tag, ".pdf"))
plot_width <- 6
plot_height <- 5
ggsave(out_file, p, width = plot_width, height = plot_height)
cat("Saved plot ->", out_file, "\n")
