rm(list = ls())
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

BASE_DIR    <- "/proj/c.zihao/work1/02analysis/04time"
RESULTS_DIR <- file.path(BASE_DIR, "results")
PLOT_DIR    <- file.path(BASE_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

py_file <- file.path(RESULTS_DIR, "timing_results_py.csv")
r_file  <- file.path(RESULTS_DIR, "timing_results_LIONESS.csv")

method_levels <- c(
  "MOSSN", "SSN", "SWEET", "Patkar", "Proteinarium", "PPIXpress", "LIONESS"
)

method_colors <- c(
  MOSSN = "#b2182b",
  SSN = "#82a1c2",
  SWEET = "#9ccdc9",
  Patkar = "#f8aa5d",
  Proteinarium = "#ba9e8e",
  PPIXpress = "#c9a1be",
  LIONESS = "#87be81"
)

read_merged <- function(merged_file) {
  if (!file.exists(merged_file)) {
    stop("Required aggregated timing table not found: ", merged_file,
         ". Run 03aggregate.R first.")
  }
  fread(merged_file)
}

py <- read_merged(py_file)
li <- read_merged(r_file)

df <- rbind(py, li, fill = TRUE)

df[, wall_time_s_chr := as.character(wall_time_s)]
df <- df[!is.na(wall_time_s_chr) & wall_time_s_chr != "ERROR"]
df[, wall_time_s       := as.numeric(wall_time_s)]
df[, time_per_sample_s := as.numeric(time_per_sample_s)]
df[, peak_memory_mb    := as.numeric(peak_memory_mb)]
df <- df[!is.na(time_per_sample_s) & time_per_sample_s > 0]

sum_df <- df[, .(
  avg_time_per_sample = mean(time_per_sample_s, na.rm = TRUE),
  avg_wall_time       = mean(wall_time_s,       na.rm = TRUE),
  avg_peak_mb         = mean(peak_memory_mb,    na.rm = TRUE)
), by = .(method, n_samples, n_edges)]

sample_df <- sum_df[n_edges == 10000L]
edge_df   <- sum_df[n_samples == 10L]

present_methods <- unique(as.character(sum_df$method))
method_order <- method_levels[method_levels %in% present_methods]
if (length(method_order) == 0) {
  method_order <- sort(present_methods)
}

sample_df[, method := factor(method, levels = method_order)]
edge_df[, method := factor(method, levels = method_order)]
sum_df[, method := factor(method, levels = method_order)]

present_colors <- method_colors[names(method_colors) %in% method_order]

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_line_samples <- ggplot(sample_df, aes(x = n_samples, y = avg_time_per_sample,
                              colour = method, group = method)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.8) +
  scale_x_log10(breaks = c(10, 20, 50, 100)) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(values = present_colors, drop = FALSE) +
  labs(x = "Number of samples", y = "Time per sample (s)",
       colour = "Method",
       title = "Computational efficiency vs sample size\n(fixed edges = 10k)") +
  base_theme

p_line_edges <- ggplot(edge_df, aes(x = n_edges, y = avg_time_per_sample,
                                    colour = method, group = method)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.8) +
  scale_x_log10(breaks = c(10000, 20000, 50000, 100000),
               labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(values = present_colors, drop = FALSE) +
  labs(x = "Number of edges", y = "Time per sample (s)",
       colour = "Method",
       title = "Computational efficiency vs network size\n(fixed samples = 10)") +
  base_theme

bar_df <- sample_df[n_samples == 100L]
bar_df <- bar_df[order(avg_time_per_sample, decreasing = TRUE)]
bar_df[, method := factor(method, levels = as.character(method))]

BAR_FLOOR <- 0.01

p_bar <- ggplot(bar_df, aes(x = method, colour = method)) +
  geom_segment(aes(xend = method, y = BAR_FLOOR, yend = avg_time_per_sample),
               linewidth = 0.5) +
  geom_point(aes(y = avg_time_per_sample), size = 1.2) +
  coord_flip() +
  scale_y_log10(limits = c(BAR_FLOOR, NA),
                labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(values = present_colors, drop = FALSE) +
  labs(x = NULL, y = "Time per sample (s)",
       title = "Time per sample at N = 100, E = 10k") +
  base_theme +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 5),
    plot.title = element_text(size = 5.5),
    plot.title.position = "plot",
    plot.margin = margin(2, 4, 2, 2)
  )

ggsave(file.path(PLOT_DIR, "figure3f.pdf"),
       p_bar, width = 5, height = 5, units = "cm")

mem_df <- sample_df[n_samples == 100L & !is.na(avg_peak_mb)]

p_mem <- NULL
if (nrow(mem_df) > 0) {
  mem_df <- mem_df[order(avg_peak_mb, decreasing = TRUE)]
  mem_df[, method := factor(method, levels = as.character(method))]

  p_mem <- ggplot(mem_df, aes(x = method, y = avg_peak_mb, fill = method)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = present_colors, drop = FALSE) +
    labs(x = NULL, y = "Peak memory (MB)",
         title = "Peak memory at N = 100, E = 10k\n(Python methods)") +
    base_theme +
    theme(legend.position = "none")
}

small_theme <- theme(
  text = element_text(size = 6),
  axis.text = element_text(size = 5),
  plot.title = element_text(size = 5.5),
  plot.title.position = "plot",
  plot.tag = element_text(size = 7, face = "bold"),
  legend.key.size = unit(0.3, "cm"),
  plot.margin = margin(2, 3, 2, 3)
)

combined <- if (!is.null(p_mem)) {
  (p_line_samples | p_line_edges | p_mem) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") & small_theme
} else {
  (p_line_samples | p_line_edges) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") & small_theme
}

ggsave(file.path(PLOT_DIR, "supfigure9.pdf"),
       combined, width = 18, height = 6, units = "cm")

sample_wide <- dcast(sample_df, method ~ n_samples,
                     value.var = "avg_time_per_sample")
existing_sample_cols <- intersect(colnames(sample_wide), as.character(c(10, 20, 50, 100)))
setnames(sample_wide, existing_sample_cols,
         paste0("N", existing_sample_cols, "_s_per_sample"))
fwrite(sample_wide, file.path(RESULTS_DIR, "timing_summary_samples.csv"))

edge_wide <- dcast(edge_df, method ~ n_edges,
                   value.var = "avg_time_per_sample")
existing_edge_cols <- intersect(colnames(edge_wide), as.character(c(10000L, 20000L, 50000L, 100000L)))
setnames(edge_wide, existing_edge_cols,
         paste0("E", existing_edge_cols, "_s_per_sample"))
fwrite(edge_wide, file.path(RESULTS_DIR, "timing_summary_edges.csv"))

cat("Plots saved to", PLOT_DIR, "\n")
print(sample_wide)
print(edge_wide)
