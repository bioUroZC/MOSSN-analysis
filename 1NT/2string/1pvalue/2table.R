rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(flextable)

out_dir <- paste0(PROJ_ROOT, "/1NT/2string/1pvalue")

df <- read.csv(file.path(out_dir, "pvalue.csv"), stringsAsFactors = FALSE)

# friendlier column names / labels for the rendered table
colnames(df) <- c("Comparison", "Group", "Method", "Accuracy p", "AUC p")
df$Comparison <- ifelse(df$Comparison == "benchmark", "vs. benchmark", "vs. ablation")
df[["Accuracy p"]] <- sprintf("%.3f", df[["Accuracy p"]])
df[["AUC p"]] <- sprintf("%.3f", df[["AUC p"]])

ft <- flextable(df)
ft <- theme_booktabs(ft)
ft <- bold(ft, part = "header")
ft <- merge_v(ft, j = c("Comparison", "Group"))
ft <- valign(ft, j = c("Comparison", "Group"), valign = "top")
ft <- align(ft, j = c("Accuracy p", "AUC p"), align = "center", part = "all")
ft <- align(ft, j = c("Comparison", "Group", "Method"), align = "left", part = "all")
ft <- font(ft, fontname = "Helvetica", part = "all")
ft <- fontsize(ft, size = 8, part = "all")
ft <- padding(ft, padding = 3, part = "all")
ft <- add_header_lines(ft, values = "MOSSN vs. benchmark methods & ablation variants (one-sided paired Wilcoxon test on cancer-type-level means, n = 11, accuracy / AUC)")
ft <- autofit(ft)

gr <- gen_grob(ft, fit = "auto")
w  <- max(11, flextable_dim(ft)$widths)
h  <- max(3, flextable_dim(ft)$heights)

cairo_pdf(file.path(out_dir, "pvalue_table.pdf"), width = w, height = h)
grid::grid.draw(gr); dev.off()

save_as_docx(ft, path = file.path(out_dir, "pvalue_table.docx"))

cat("Rendered pvalue table -> pdf / docx\n")
