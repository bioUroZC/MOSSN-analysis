rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
module_dir <- file.path(base_dir, "module_results")

input_file <- file.path(module_dir, "module_summary.csv")
out_pdf <- file.path(base_dir, "module_gene_counts.pdf")
out_csv <- file.path(base_dir, "module_gene_counts.csv")

module_summary <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

plot_tbl <- module_summary %>%
  dplyr::arrange(direction, dplyr::desc(n_genes), module_id) %>%
  dplyr::mutate(
    module_id = factor(module_id, levels = rev(module_id)),
    direction = factor(direction, levels = c("gain", "loss"))
  )

utils::write.csv(plot_tbl, out_csv, row.names = FALSE)

p <- ggplot(plot_tbl, aes(x = n_genes, y = module_id, fill = direction)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n_genes), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Number of genes",
    y = "Module",
    fill = "Direction",
    title = "Gene counts across recurrent rewiring modules"
  ) +
  expand_limits(x = max(plot_tbl$n_genes, na.rm = TRUE) * 1.12)

pdf(out_pdf, width = 8.5, height = 10)
print(p)
dev.off()

message("Saved: ", out_pdf)
message("Saved: ", out_csv)
