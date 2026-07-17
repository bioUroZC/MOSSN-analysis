rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

BASE_DIR    <- paste0(PROJ_ROOT, "/1NT/6time")
RESULTS_DIR <- file.path(BASE_DIR, "results")
PLOT_DIR    <- file.path(BASE_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

py_file <- file.path(RESULTS_DIR, "timing_results_py.csv")
r_file  <- file.path(RESULTS_DIR, "timing_results_LIONESS.csv")

method_levels <- c(
  "MOSS_noPrior", "MOSSN_uniform", "MOSS_uniform", "MOSS_noSeed", "MOSS_noRWR",
  "MOSS_noCorr", "SSN", "SWEET", "Patkar", "Proteinarium", "PPIXpress", "LIONESS"
)

method_colors <- c(
  MOSS_noPrior = "#8c1c13",
  MOSSN_uniform = "#b22222",
  MOSS_uniform = "#ef8a62",
  MOSS_noSeed = "#f4a259",
  MOSS_noRWR = "#ffd166",
  MOSS_noCorr = "#3d405b",
  SSN = "#4c78a8",
  SWEET = "#72b7b2",
  Patkar = "#e4572e",
  Proteinarium = "#9d755d",
  PPIXpress = "#b279a2",
  LIONESS = "#54a24b"
)

read_or_combine <- function(merged_file, pattern) {
  if (file.exists(merged_file)) {
    return(fread(merged_file))
  }
  files <- Sys.glob(file.path(RESULTS_DIR, pattern))
  if (length(files) == 0) {
    stop("No timing result files found for pattern: ", pattern)
  }
  rbindlist(lapply(sort(files), fread), fill = TRUE)
}

# ── 1. Load and merge Python + LIONESS results ────────────────────────────────
py <- read_or_combine(py_file, "timing_S*_E*.csv")
li <- read_or_combine(r_file, "timing_LIONESS_S*_E*.csv")

df <- rbind(py, li, fill = TRUE)

# Keep only successful numeric rows
df[, wall_time_s_chr := as.character(wall_time_s)]
df <- df[!is.na(wall_time_s_chr) & wall_time_s_chr != "ERROR"]
df[, wall_time_s       := as.numeric(wall_time_s)]
df[, time_per_sample_s := as.numeric(time_per_sample_s)]
df[, peak_memory_mb    := as.numeric(peak_memory_mb)]
df <- df[!is.na(time_per_sample_s) & time_per_sample_s > 0]

# ── 2. Summarise: median over reps per grid point ─────────────────────────────
sum_df <- df[, .(
  med_time_per_sample = median(time_per_sample_s, na.rm = TRUE),
  med_wall_time       = median(wall_time_s,       na.rm = TRUE),
  med_peak_mb         = median(peak_memory_mb,    na.rm = TRUE),
  q25_time            = quantile(time_per_sample_s, 0.25, na.rm = TRUE),
  q75_time            = quantile(time_per_sample_s, 0.75, na.rm = TRUE)
), by = .(method, n_samples, n_edges)]

sample_df <- sum_df[n_edges == 10000L]
edge_df   <- sum_df[n_samples == 10L]

# Keep the same method ordering and palette used in the other 1NT figures.
present_methods <- unique(as.character(sum_df$method))
method_order <- method_levels[method_levels %in% present_methods]
if (length(method_order) == 0) {
  method_order <- sort(present_methods)
}

# MOSS vs comparison colour group
moss_methods <- grep("^MOSS", method_order, value = TRUE)
sample_df[, group := ifelse(method %in% moss_methods, "MOSS variant", "Comparison")]
edge_df[, group := ifelse(method %in% moss_methods, "MOSS variant", "Comparison")]
sum_df[, group := ifelse(method %in% moss_methods, "MOSS variant", "Comparison")]

sample_df[, method := factor(method, levels = method_order)]
edge_df[, method := factor(method, levels = method_order)]
sum_df[, method := factor(method, levels = method_order)]

present_colors <- method_colors[names(method_colors) %in% method_order]

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ── 3. Plot A: time per sample vs N at fixed edge count ───────────────────────
p_line_samples <- ggplot(sample_df, aes(x = n_samples, y = med_time_per_sample,
                              colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = q25_time, ymax = q75_time, fill = method),
              alpha = 0.10, colour = NA) +
  scale_x_log10(breaks = c(10, 20, 50, 100)) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(values = present_colors, drop = FALSE) +
  scale_fill_manual(values = alpha(present_colors, 0.18), drop = FALSE) +
  labs(x = "Number of samples", y = "Time per sample (s)",
       colour = "Method", fill = "Method",
       title = "Computational efficiency vs sample size (fixed edges = 10k)") +
  base_theme

# ── 4. Plot B: time per sample vs edge count at fixed sample size ─────────────
p_line_edges <- ggplot(edge_df, aes(x = n_edges, y = med_time_per_sample,
                                    colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = q25_time, ymax = q75_time, fill = method),
              alpha = 0.10, colour = NA) +
  scale_x_log10(breaks = c(10000, 20000, 50000, 100000),
               labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(values = present_colors, drop = FALSE) +
  scale_fill_manual(values = alpha(present_colors, 0.18), drop = FALSE) +
  labs(x = "Number of edges", y = "Time per sample (s)",
       colour = "Method", fill = "Method",
       title = "Computational efficiency vs network size (fixed samples = 10)") +
  base_theme

# ── 5. Plot C: bar chart at N=100, E=10k ─────────────────────────────────────
bar_df <- sample_df[n_samples == 100L]
bar_df <- bar_df[order(med_time_per_sample, decreasing = TRUE)]
bar_df[, method := factor(method, levels = as.character(method))]

p_bar <- ggplot(bar_df, aes(x = method, y = med_time_per_sample, fill = method)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = q25_time, ymax = q75_time), width = 0.25) +
  coord_flip() +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  scale_fill_manual(values = present_colors, drop = FALSE) +
  labs(x = NULL, y = "Time per sample (s)",
       title = "Time per sample at N = 100, E = 10k") +
  base_theme +
  theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "timing_bar_N100.pdf"),
       p_bar, width = 6, height = 5)

# ── 6. Plot D: peak memory bar at N=100, E=10k (Python methods only) ─────────
mem_df <- sample_df[n_samples == 100L & !is.na(med_peak_mb)]

p_mem <- NULL
if (nrow(mem_df) > 0) {
  mem_df <- mem_df[order(med_peak_mb, decreasing = TRUE)]
  mem_df[, method := factor(method, levels = as.character(method))]

  p_mem <- ggplot(mem_df, aes(x = method, y = med_peak_mb, fill = method)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = present_colors, drop = FALSE) +
    labs(x = NULL, y = "Peak memory (MB)",
         title = "Peak memory at N = 100, E = 10k (Python methods)") +
    base_theme +
    theme(legend.position = "none")
}

# ── 6b. Combined PDF: three plots in one row ──────────────────────────────────
combined <- if (!is.null(p_mem)) {
  (p_line_samples | p_line_edges | p_mem) +
    plot_annotation(tag_levels = "A")
} else {
  (p_line_samples | p_line_edges) +
    plot_annotation(tag_levels = "A")
}

ggsave(file.path(PLOT_DIR, "timing_combined.pdf"),
       combined, width = 18, height = 5)

# ── 7. Save summary tables ────────────────────────────────────────────────────
sample_wide <- dcast(sample_df, method ~ n_samples,
                     value.var = "med_time_per_sample")
existing_sample_cols <- intersect(colnames(sample_wide), as.character(c(10, 20, 50, 100)))
setnames(sample_wide, existing_sample_cols,
         paste0("N", existing_sample_cols, "_s_per_sample"))
fwrite(sample_wide, file.path(RESULTS_DIR, "timing_summary_samples.csv"))

edge_wide <- dcast(edge_df, method ~ n_edges,
                   value.var = "med_time_per_sample")
existing_edge_cols <- intersect(colnames(edge_wide), as.character(c(10000, 20000, 50000, 100000)))
setnames(edge_wide, existing_edge_cols,
         paste0("E", existing_edge_cols, "_s_per_sample"))
fwrite(edge_wide, file.path(RESULTS_DIR, "timing_summary_edges.csv"))

cat("Plots saved to", PLOT_DIR, "\n")
print(sample_wide)
print(edge_wide)
