rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

base_family <- "Helvetica"

IN_DIR  <- paste0(PROJ_ROOT, "/4coupled/noise/results/spearman_noise")
OUT_DIR <- file.path(IN_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

METHODS <- c("MOSSN_EXP", "MOSSN_Direct", "MOSSN_DirectNoDyn", "MOSSN_Restart", "MOSSN_MultiLayer")
NOISE_LEVELS <- c("k0.0", "k0.5", "k1.0", "k1.5", "k2.0", "k3.0", "k5.0")
NOISE_LABELS <- c("0", "0.5", "1.0", "1.5", "2.0", "3.0", "5.0")

METHOD_LABELS <- c(
    "MOSSN_EXP" = "MOSSN-EXP",
    "MOSSN_Direct" = "MOSSN-Direct",
    "MOSSN_DirectNoDyn" = "MOSSN-DirectNoDyn",
    "MOSSN_Restart" = "MOSSN-Restart",
    "MOSSN_MultiLayer" = "MOSSN-MultiLayer"
)

METHOD_COLS <- c(
    "MOSSN_EXP"            = "#1b9e77",
    "MOSSN_Direct"         = "#d7301f",
    "MOSSN_DirectNoDyn"    = "#ef6548",
    "MOSSN_Restart"        = "#2c7fb8",
    "MOSSN_MultiLayer"     = "#c51b7d"
)

NOISE_NUM <- c("k0.5" = 0.5, "k1.0" = 1.0, "k1.5" = 1.5, "k2.0" = 2.0, "k3.0" = 3.0, "k5.0" = 5.0)

METHOD_COLS_LABELED <- setNames(
    unname(METHOD_COLS[METHODS]),
    METHOD_LABELS[METHODS]
)

base_theme <- theme_bw(base_size = 12, base_family = base_family) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
        panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6),
        axis.title = element_text(face = "bold", family = base_family),
        axis.text = element_text(color = "#222222", family = base_family),
        plot.title = element_text(face = "bold", size = 14, family = base_family),
        plot.subtitle = element_text(size = 10.5, color = "#4d4d4d", family = base_family),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, family = base_family),
        legend.spacing.x = grid::unit(4, "pt"),
        legend.key.width = grid::unit(14, "pt"),
        plot.margin = margin(8, 10, 8, 8)
    )

df <- read.csv(file.path(IN_DIR, "spearman_results.csv"), stringsAsFactors = FALSE)
df <- df[df$Method %in% METHODS, ]
df$Method    <- factor(df$Method, levels = METHODS)
df$NoiseNum  <- NOISE_NUM[df$NoiseLevel]

# add k=0 baseline (rho=1 by definition)
baseline <- expand.grid(Method = METHODS, NoiseLevel = "k0.0", stringsAsFactors = FALSE)
baseline$NoiseNum       <- 0
baseline$Median_Spearman <- 1
baseline$Q25_Spearman   <- 1
baseline$Q75_Spearman   <- 1
baseline$Method <- factor(baseline$Method, levels = METHODS)

plot_df <- bind_rows(df, baseline)
plot_df$Method <- factor(plot_df$Method, levels = METHODS)
plot_df$NoiseLevel <- factor(plot_df$NoiseLevel, levels = NOISE_LEVELS)
plot_df$MethodLabel <- factor(
    METHOD_LABELS[as.character(plot_df$Method)],
    levels = unname(METHOD_LABELS[METHODS])
)
plot_df$Highlight <- ifelse(plot_df$Method == "MOSSN_Direct", "highlight", "background")

# line + ribbon (IQR)
p <- ggplot() +
    geom_ribbon(
        data = plot_df[plot_df$Highlight == "background", ],
        aes(
            x = NoiseLevel,
            ymin = Q25_Spearman,
            ymax = Q75_Spearman,
            fill = MethodLabel,
            group = MethodLabel
        ),
        alpha = 0.10,
        color = NA
    ) +
    geom_ribbon(
        data = plot_df[plot_df$Highlight == "highlight", ],
        aes(
            x = NoiseLevel,
            ymin = Q25_Spearman,
            ymax = Q75_Spearman,
            fill = MethodLabel,
            group = MethodLabel
        ),
        alpha = 0.14,
        color = NA
    ) +
    geom_line(
        data = plot_df[plot_df$Highlight == "background", ],
        aes(x = NoiseLevel, y = Median_Spearman, color = MethodLabel, group = MethodLabel),
        linewidth = 0.9,
        alpha = 0.95
    ) +
    geom_point(
        data = plot_df[plot_df$Highlight == "background", ],
        aes(x = NoiseLevel, y = Median_Spearman, color = MethodLabel),
        size = 2.6,
        alpha = 0.98
    ) +
    geom_line(
        data = plot_df[plot_df$Highlight == "highlight", ],
        aes(x = NoiseLevel, y = Median_Spearman, color = MethodLabel, group = MethodLabel),
        linewidth = 1.45
    ) +
    geom_point(
        data = plot_df[plot_df$Highlight == "highlight", ],
        aes(x = NoiseLevel, y = Median_Spearman, color = MethodLabel),
        size = 3.1
    ) +
    scale_x_discrete(labels = NOISE_LABELS, drop = FALSE) +
    scale_y_continuous(limits = c(NA, 1), expand = expansion(mult = c(0.02, 0.03))) +
    scale_color_manual(values = METHOD_COLS_LABELED, drop = FALSE) +
    scale_fill_manual(values = METHOD_COLS_LABELED, drop = FALSE) +
    base_theme +
    labs(
        title    = "Noise robustness of network consistency in LUAD",
        subtitle = "Median Spearman correlation between noisy-edge and original-edge networks; ribbons indicate IQR",
        x        = "Noise level k",
        y        = "Median Spearman rho"
    )

ggsave(file.path(OUT_DIR, "noise_spearman.pdf"), p, width = 8.6, height = 5.8)
cat("Saved ->", file.path(OUT_DIR, "noise_spearman.pdf"), "\n")

# delta vs MOSSN_EXP at each non-zero noise level
delta_df <- df |>
    select(Method, NoiseLevel, NoiseNum, Median_Spearman, Q25_Spearman, Q75_Spearman) |>
    left_join(
        df |>
            filter(Method == "MOSSN_EXP") |>
            select(
                NoiseLevel,
                Base_Median = Median_Spearman,
                Base_Q25 = Q25_Spearman,
                Base_Q75 = Q75_Spearman
            ),
        by = "NoiseLevel"
    ) |>
    filter(Method != "MOSSN_EXP") |>
    mutate(
        Delta_Median = Median_Spearman - Base_Median,
        Delta_Q25 = Q25_Spearman - Base_Q25,
        Delta_Q75 = Q75_Spearman - Base_Q75
    )

delta_df$NoiseLevel <- factor(delta_df$NoiseLevel, levels = NOISE_LEVELS[-1])
delta_df$Method <- factor(delta_df$Method, levels = METHODS[METHODS != "MOSSN_EXP"])
delta_df$MethodLabel <- factor(
    METHOD_LABELS[as.character(delta_df$Method)],
    levels = unname(METHOD_LABELS[METHODS[METHODS != "MOSSN_EXP"]])
)
delta_df$Highlight <- ifelse(delta_df$Method == "MOSSN_Direct", "highlight", "background")

p_delta <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#8c8c8c", linewidth = 0.5) +
    geom_line(
        data = delta_df[delta_df$Highlight == "background", ],
        aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel, group = MethodLabel),
        linewidth = 0.9,
        alpha = 0.95
    ) +
    geom_point(
        data = delta_df[delta_df$Highlight == "background", ],
        aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel),
        size = 2.6,
        alpha = 0.98
    ) +
    geom_line(
        data = delta_df[delta_df$Highlight == "highlight", ],
        aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel, group = MethodLabel),
        linewidth = 1.45
    ) +
    geom_point(
        data = delta_df[delta_df$Highlight == "highlight", ],
        aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel),
        size = 3.1
    ) +
    scale_x_discrete(labels = NOISE_LABELS[-1], drop = FALSE) +
    scale_color_manual(
        values = METHOD_COLS_LABELED[unname(METHOD_LABELS[METHODS[METHODS != "MOSSN_EXP"]])],
        drop = FALSE
    ) +
    base_theme +
    labs(
        x        = "Noise level k",
        y        = "Delta Spearman rho"
    )

ggsave(file.path(OUT_DIR, "noise_spearman_delta.pdf"), p_delta, width = 6.2, height = 6.2)
cat("Saved ->", file.path(OUT_DIR, "noise_spearman_delta.pdf"), "\n")
