rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

BASE_DIR <- "/proj/c.zihao/work1/04multi"
OUT_FILE <- file.path(BASE_DIR, "supfigure11.pdf")

NOISE_LEVELS <- c("k0.5", "k1.0", "k1.5", "k2.0", "k3.0", "k5.0")
NOISE_LABELS <- c("0.5", "1.0", "1.5", "2.0", "3.0", "5.0")

METHOD_LABELS <- c(
  Direct = "Direct",
  NoDyn = "NoDyn",
  Restart = "Restart",
  MultiLayer = "MultiLayer"
)

METHOD_COLS <- c(
  Direct = "#d7301f",
  NoDyn = "#ef6548",
  Restart = "#2c7fb8",
  MultiLayer = "#c51b7d"
)

base_theme <- theme_bw(base_size = 3, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
    panel.border = element_rect(color = "#4d4d4d", linewidth = 0.6),
    axis.title = element_text(size = 3, face = "bold", family = "Helvetica"),
    axis.text = element_text(size = 3, color = "#222222", family = "Helvetica"),
    plot.subtitle = element_text(size = 3, family = "Helvetica"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 3, family = "Helvetica"),
    legend.spacing.x = grid::unit(2, "pt"),
    legend.key.width = grid::unit(8, "pt"),
    legend.key.height = grid::unit(3, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0, "pt"),
    plot.margin = margin(3, 4, 3, 3)
  )

make_delta_plot <- function(result_dir, baseline_method, subtitle) {
  methods <- c(baseline_method, names(METHOD_LABELS))

  df <- read.csv(
    file.path(result_dir, "spearman.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df <- df[df$Method %in% methods & df$NoiseLevel %in% NOISE_LEVELS, ]

  delta_df <- df %>%
    select(Method, NoiseLevel, Median_Spearman, Q25_Spearman, Q75_Spearman) %>%
    left_join(
      df %>%
        filter(Method == baseline_method) %>%
        select(
          NoiseLevel,
          Base_Median = Median_Spearman,
          Base_Q25 = Q25_Spearman,
          Base_Q75 = Q75_Spearman
        ),
      by = "NoiseLevel"
    ) %>%
    filter(Method != baseline_method) %>%
    mutate(
      Delta_Median = Median_Spearman - Base_Median,
      Delta_Q25 = Q25_Spearman - Base_Q25,
      Delta_Q75 = Q75_Spearman - Base_Q75,
      MethodLabel = factor(Method, levels = names(METHOD_LABELS)),
      Highlight = ifelse(Method == "Direct", "highlight", "background")
    )

  ggplot() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "#8c8c8c",
      linewidth = 0.2
    ) +
    geom_line(
      data = delta_df[delta_df$Highlight == "background", ],
      aes(
        x = NoiseLevel,
        y = Delta_Median,
        color = MethodLabel,
        group = MethodLabel
      ),
      linewidth = 0.3,
      alpha = 0.95
    ) +
    geom_point(
      data = delta_df[delta_df$Highlight == "background", ],
      aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel),
      size = 0.9,
      alpha = 0.98
    ) +
    geom_line(
      data = delta_df[delta_df$Highlight == "highlight", ],
      aes(
        x = NoiseLevel,
        y = Delta_Median,
        color = MethodLabel,
        group = MethodLabel
      ),
      linewidth = 0.5
    ) +
    geom_point(
      data = delta_df[delta_df$Highlight == "highlight", ],
      aes(x = NoiseLevel, y = Delta_Median, color = MethodLabel),
      size = 1.2
    ) +
    scale_x_discrete(
      limits = NOISE_LEVELS,
      labels = NOISE_LABELS,
      drop = FALSE
    ) +
    scale_color_manual(
      values = METHOD_COLS,
      labels = METHOD_LABELS,
      drop = FALSE
    ) +
    base_theme +
    labs(
      subtitle = subtitle,
      x = "Noise level k",
      y = "Delta Spearman rho"
    )
}

p_cnv <- make_delta_plot(
  file.path(BASE_DIR, "3cnvNoise/04result"),
  "CNV",
  "Robustness of MOSSN ablation variants\nunder increasing CNV noise"
)
p_exp <- make_delta_plot(
  file.path(BASE_DIR, "3expNoise/04result"),
  "EXP",
  "Robustness of MOSSN ablation variants\nunder increasing expression noise"
)
p_met <- make_delta_plot(
  file.path(BASE_DIR, "3metNoise/04result"),
  "MET",
  "Robustness of MOSSN ablation variants\nunder increasing methylation noise"
)
p_all <- make_delta_plot(
  file.path(BASE_DIR, "4allNoise/04result"),
  "EXP",
  "Robustness of MOSSN ablation variants\nunder increasing all-layer noise"
)

combined <- ((p_cnv + p_exp) / (p_met + p_all)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 3, face = "bold"))
  )
combined <- combined & theme(legend.position = "bottom")

ggsave(
  filename = OUT_FILE,
  plot = combined,
  width = 8,
  height = 8,
  units = "cm"
)

cat("Saved ->", OUT_FILE, "\n")
