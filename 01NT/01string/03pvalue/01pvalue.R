rm(list = ls())

base_dir <- "/proj/c.zihao/work1/01NT/01string"
benchmark_csv <- file.path(base_dir, "02benchmark/02merge/results_Cluster.csv")
ablation_csv <- file.path(base_dir, "01ablation/02merge/results_Cluster.csv")
output_dir <- file.path(base_dir, "03pvalue")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

full_method <- "MOSSN"
metrics <- c("accuracy", "auc")

benchmark_method_levels <- c(
  "SSN", "SWEET", "LIONESS", "Patkar", "Proteinarium", "PPIXpress"
)
ablation_method_levels <- c(
  "noCorr", "noRWR", "noSeed",
  "RandomBackbone", "PermutedControl",
  "NodeRWR", "RawExpr"
)

paired_compare <- function(df, fractions = NULL, group_label = "all") {
  if (!is.null(fractions)) {
    df <- df[df$feature_fraction %in% fractions, ]
  }

  agg <- aggregate(
    cbind(accuracy, auc) ~ cancer + method,
    data = df,
    FUN = mean
  )

  full_df <- agg[agg$method == full_method, ]
  other_methods <- setdiff(unique(agg$method), full_method)

  out <- data.frame(
    group = character(),
    method = character(),
    metric = character(),
    n = integer(),
    mossn_mean = numeric(),
    other_mean = numeric(),
    mean_diff = numeric(),
    wilcox_p = numeric(),
    stringsAsFactors = FALSE
  )

  for (m in other_methods) {
    other_df <- agg[agg$method == m, ]

    for (metric in metrics) {
      merged <- merge(
        full_df[, c("cancer", metric)],
        other_df[, c("cancer", metric)],
        by = "cancer",
        suffixes = c(".mossn", ".other")
      )
      mossn_vals <- merged[[paste0(metric, ".mossn")]]
      other_vals <- merged[[paste0(metric, ".other")]]
      n <- nrow(merged)

      wilcox_p <- tryCatch(
        wilcox.test(mossn_vals, other_vals, paired = TRUE,
                     alternative = "greater")$p.value,
        error = function(e) NA_real_
      )

      out <- rbind(out, data.frame(
        group = group_label,
        method = m,
        metric = metric,
        n = n,
        mossn_mean = round(mean(mossn_vals), 5),
        other_mean = round(mean(other_vals), 5),
        mean_diff = round(mean(mossn_vals - other_vals), 5),
        wilcox_p = signif(wilcox_p, 4),
        stringsAsFactors = FALSE
      ))
    }
  }

  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_real_, ifelse(p < 0.001, 0.001, round(p, 3)))
}

benchmark_df <- read.csv(benchmark_csv, stringsAsFactors = FALSE)
benchmark_result <- paired_compare(benchmark_df, fractions = NULL, group_label = "all")
benchmark_result <- cbind(comparison = "benchmark", benchmark_result)

ablation_df <- read.csv(ablation_csv, stringsAsFactors = FALSE)

mossn_rows <- benchmark_df[benchmark_df$method == full_method, ]
ablation_full_df <- rbind(ablation_df, mossn_rows[, colnames(ablation_df)])

ablation_high <- paired_compare(ablation_full_df, fractions = c("15%", "20%"), group_label = "high(15-20%)")
ablation_low <- paired_compare(ablation_full_df, fractions = c("5%", "10%"), group_label = "low(5-10%)")

ablation_result <- rbind(ablation_high, ablation_low)
ablation_result <- cbind(comparison = "ablation", ablation_result)

pvalue_long <- rbind(benchmark_result, ablation_result)

wide <- reshape(
  pvalue_long[, c("comparison", "group", "method", "metric", "n", "wilcox_p")],
  idvar = c("comparison", "group", "method"),
  timevar = "metric",
  direction = "wide"
)

pvalue_result <- data.frame(
  comparison = wide$comparison,
  group = wide$group,
  method = wide$method,
  accuracy_p = fmt_p(wide$wilcox_p.accuracy),
  auc_p = fmt_p(wide$wilcox_p.auc),
  stringsAsFactors = FALSE
)
method_order <- factor(
  pvalue_result$method,
  levels = c(benchmark_method_levels, ablation_method_levels)
)
comparison_order <- factor(pvalue_result$comparison, levels = c("benchmark", "ablation"))
group_order <- factor(pvalue_result$group, levels = c("all", "low(5-10%)", "high(15-20%)"))
pvalue_result <- pvalue_result[order(comparison_order, group_order, method_order), ]

pvalue_result$comparison <- ifelse(
  pvalue_result$comparison == "benchmark", "vs. benchmark", "vs. ablation"
)
pvalue_result$group <- ifelse(
  pvalue_result$group == "low(5-10%)", "Low (5-10%)",
  ifelse(pvalue_result$group == "high(15-20%)", "High(15-20%)", pvalue_result$group)
)
colnames(pvalue_result) <- c(
  "Comparison", "Group", "Method", "P-value (Accuracy)", "P-value (AUC)"
)

cat("\n=== MOSSN vs benchmark methods & ablation variants (Wilcoxon paired p) ===\n")
print(pvalue_result)

write.csv(
  pvalue_result,
  file.path(output_dir, "Supplementary Table 5.csv"),
  row.names = FALSE
)

cat("\nSaved:", file.path(output_dir, "Supplementary Table 5.csv"), "\n")
