rm(list = ls())
library(flextable)

out_dir <- "/proj/c.zihao/work1/01NT/01string/03pvalue"

df <- read.csv(
  file.path(out_dir, "Supplementary Table 5.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df$Group <- ifelse(
  df$Group == "low(5-10%)", "Low (5-10%)",
  ifelse(df$Group == "high(15-20%)", "High(15-20%)", df$Group)
)
df[["P-value (Accuracy)"]] <- sprintf("%.3f", df[["P-value (Accuracy)"]])
df[["P-value (AUC)"]] <- sprintf("%.3f", df[["P-value (AUC)"]])

ft <- flextable(df)
ft <- theme_booktabs(ft)
ft <- bold(ft, part = "header")
ft <- merge_v(ft, j = c("Comparison", "Group"))
ft <- valign(ft, j = c("Comparison", "Group"), valign = "top")
ft <- align(ft, j = c("P-value (Accuracy)", "P-value (AUC)"), align = "center", part = "all")
ft <- align(ft, j = c("Comparison", "Group", "Method"), align = "left", part = "all")
ft <- font(ft, fontname = "Helvetica", part = "all")
ft <- fontsize(ft, size = 8, part = "all")
ft <- padding(ft, padding = 3, part = "all")
ft <- add_header_lines(ft, values = "MOSSN vs. benchmark methods & ablation variants (one-sided paired Wilcoxon test on cancer-type-level means, n = 11, accuracy / AUC)")
ft <- autofit(ft)

save_as_docx(ft, path = file.path(out_dir, "Supplementary Table 5.docx"))

cat("Rendered Supplementary Table 5 -> docx\n")
