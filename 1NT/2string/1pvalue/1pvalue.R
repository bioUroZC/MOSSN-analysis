rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


base_dir <- paste0(PROJ_ROOT, "/1NT/2string")
benchmark_csv <- file.path(base_dir, "benchmark/1result/results_Cluster.csv")
ablation_csv <- file.path(base_dir, "ablation/1result/results_Cluster.csv")
output_dir <- file.path(base_dir, "1pvalue")

full_method <- "MOSSN_uniform"
metrics <- c("accuracy", "auc")

benchmark_method_levels <- c(
  "SSN", "SWEET", "LIONESS", "Patkar", "Proteinarium", "PPIXpress"
)
ablation_method_levels <- c(
  "MOSSN_noCorr", "EdgeNoRWR", "MOSSN_noSeed",
  "RandomBackbone", "PermutedControl",
  "NodeRWR", "RawExpr"
)

# paired comparison of `full_method` against every other method present in
# `df`, at the cancer-type level. `fractions = NULL` uses all feature
# fractions; otherwise restricts to the given feature_fraction labels first.
# feature_fraction is a hyperparameter of the same patient cohort, not an
# independent replicate, so pairing on (cancer, feature_fraction) would be
# pseudoreplication (same issue already handled correctly for the survival
# analysis elsewhere in the manuscript). To avoid this, accuracy/auc are
# first averaged across the selected fractions within each cancer type, and
# the paired Wilcoxon test is run on those per-cancer means (n = #cancers).
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

      # one-sided: MOSSN is hypothesized to be no worse than the comparator
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

# round to 3 decimals; anything below 0.001 is floored to 0.001
fmt_p <- function(p) {
  ifelse(is.na(p), NA_real_, ifelse(p < 0.001, 0.001, round(p, 3)))
}

## ---- 1. MOSSN vs benchmark methods (all feature fractions) ----------------

benchmark_df <- read.csv(benchmark_csv, stringsAsFactors = FALSE)
benchmark_result <- paired_compare(benchmark_df, fractions = NULL, group_label = "all")
benchmark_result <- cbind(comparison = "benchmark", benchmark_result)

## ---- 2. MOSSN vs ablation variants (full model taken from benchmark) ------

ablation_df <- read.csv(ablation_csv, stringsAsFactors = FALSE)
ablation_df <- ablation_df[ablation_df$method != "DegreePreservedBackbone", ]

# the full MOSSN_uniform run lives in the benchmark results; splice it into
# the ablation table so every variant is compared against the same baseline.
mossn_rows <- benchmark_df[benchmark_df$method == full_method, ]
ablation_full_df <- rbind(ablation_df, mossn_rows[, colnames(ablation_df)])

ablation_high <- paired_compare(ablation_full_df, fractions = c("15%", "20%"), group_label = "high(15-20%)")
ablation_low <- paired_compare(ablation_full_df, fractions = c("5%", "10%"), group_label = "low(5-10%)")

ablation_result <- rbind(ablation_high, ablation_low)
ablation_result <- cbind(comparison = "ablation", ablation_result)

## ---- 3. merge into a single table, one row per (comparison, group, method) -

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

cat("\n=== MOSSN vs benchmark methods & ablation variants (Wilcoxon paired p) ===\n")
print(pvalue_result)

write.csv(
  pvalue_result,
  file.path(output_dir, "pvalue.csv"),
  row.names = FALSE
)

cat("\nSaved:", file.path(output_dir, "pvalue.csv"), "\n")
