rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)
library(ggplot2)
library(poweRlaw)

network  <- "STRING"
data_dir <- file.path(paste0(PROJ_ROOT, "/1NT/2free"), network)
out_dir  <- file.path(data_dir, "plots")
dir.create(out_dir, showWarnings = FALSE)

dat <- read_csv(file.path(data_dir, "scalefree_persample.csv"),
                show_col_types = FALSE)

top_levels <- c("5%", "10%", "15%", "20%")

method_order <- c(
  "Backbone (STRING)",
  "MOSSN_uniform",
  "SSN", "SWEET", "LIONESS", "Patkar",  "Proteinarium", "PPIXpress",
  "MOSSN_noCorr",  "EdgeNoRWR",  "MOSSN_noSeed",
  "RandomBackbone", "PermutedControl"
)

dat <- dat %>% filter(method %in% method_order)

# ── backbone stats per threshold ──────────────────────────────────────────────
bb_links <- read_csv(
  paste0(PROJ_ROOT, "/1NT/1data/string/links.csv"),
  show_col_types = FALSE
) %>% select(protein1, protein2, score)

bb_row <- bind_rows(lapply(c(0.05, 0.10, 0.15, 0.20), function(pct) {
  n   <- nrow(bb_links)
  k   <- max(1L, round(pct * n))
  ord <- order(bb_links$score, decreasing = TRUE)[seq_len(k)]
  deg <- as.integer(table(c(bb_links$protein1[ord], bb_links$protein2[ord])))
  deg <- deg[deg > 0L]

  pl <- tryCatch({
    o <- displ$new(deg)
    o$setXmin(estimate_xmin(o))
    o$setPars(estimate_pars(o))
    o
  }, error = function(e) NULL)

  data.frame(
    method     = "Backbone (STRING)",
    top_pct    = paste0(pct * 100, "%"),
    mean_gamma = round(if (!is.null(pl)) pl$pars else NA_real_, 3),
    mean_H     = round(mean(deg^2) / mean(deg)^2, 2),
    stringsAsFactors = FALSE
  )
}))

# ── per-method summaries ──────────────────────────────────────────────────────
sum_gamma <- dat %>%
  mutate(top_pct = paste0(top_pct * 100, "%")) %>%
  group_by(method, top_pct) %>%
  summarise(mean_gamma = mean(gamma_pl, na.rm = TRUE), .groups = "drop") %>%
  bind_rows(bb_row %>% select(method, top_pct, mean_gamma)) %>%
  mutate(
    top_pct      = factor(top_pct, levels = top_levels),
    method       = factor(method, levels = method_order),
    gamma_capped = pmin(mean_gamma, 4),
    label        = sprintf("%.2f", mean_gamma)
  )

sum_H <- dat %>%
  mutate(top_pct = paste0(top_pct * 100, "%")) %>%
  group_by(method, top_pct) %>%
  summarise(mean_H = mean(H, na.rm = TRUE), .groups = "drop") %>%
  bind_rows(bb_row %>% select(method, top_pct, mean_H)) %>%
  mutate(
    top_pct = factor(top_pct, levels = top_levels),
    method  = factor(method, levels = method_order),
    label   = sprintf("%.2f", mean_H)
  )

# ── shared theme ──────────────────────────────────────────────────────────────
theme_heatmap <- theme_bw(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    axis.text.x    = element_text(size = 9, angle = 35, hjust = 1),
    axis.text.y    = element_text(size = 10),
    legend.title   = element_text(size = 9),
    plot.title     = element_blank(),
    plot.subtitle  = element_blank()
  )

# separator between Backbone and methods (after 1st column)
vline_pos <- 1.5

# ── gamma heatmap ─────────────────────────────────────────────────────────────
p_gamma <- ggplot(sum_gamma, aes(method, top_pct, fill = gamma_capped)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.2, color = "grey10") +
  geom_vline(xintercept = vline_pos,
             color = "grey30", linewidth = 0.8) +
  scale_fill_gradientn(
    colours = c(
      "#4DAC26", "#B8E186", "#2166AC", "#2166AC", "#F4A582", "#D73027"
    ),
    values  = scales::rescale(c(1.4, 1.99, 2.0, 3.0, 3.01, 4.0)),
    limits  = c(1.4, 4),
    oob     = scales::squish,
    name    = expression(paste("Mean ", gamma)),
    guide   = guide_colorbar(barwidth = 1, barheight = 6)
  ) +
  labs(x = NULL, y = "Top-edge threshold") +
  theme_heatmap

ggsave(file.path(out_dir, "gamma_heatmap.pdf"),
       p_gamma, width = 8.6, height = 3.5)
cat("Saved: gamma_heatmap.pdf\n")

# ── H heatmap ─────────────────────────────────────────────────────────────────
p_H <- ggplot(sum_H, aes(method, top_pct, fill = mean_H)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.2, color = "grey10") +
  geom_vline(xintercept = vline_pos,
             color = "grey30", linewidth = 0.8) +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#FDD0A2", "#D94801"),
    values  = scales::rescale(c(1, 4, 10)),
    limits  = c(1, 10),
    oob     = scales::squish,
    name    = "Mean H",
    guide   = guide_colorbar(barwidth = 1, barheight = 6)
  ) +
  labs(x = NULL, y = "Top-edge threshold") +
  theme_heatmap

ggsave(file.path(out_dir, "H_heatmap.pdf"),
       p_H, width = 8.6, height = 3.5)
cat("Saved: H_heatmap.pdf\n")
