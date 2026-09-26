rm(list = ls())
suppressPackageStartupMessages(library(dplyr))

result_dir <- "/proj/c.zihao/work1/01NT/02free/result"

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

persample <- read.csv(
  file.path(result_dir, "persample.csv"),
  stringsAsFactors = FALSE
)

method_metrics <- persample %>%
  group_by(method, top_pct) %>%
  summarise(
    R2 = safe_mean(R2),
    .groups = "drop"
  ) %>%
  transmute(
    Method = method,
    top_pct,
    R2
  )

network_metrics <- function(node1, node2) {
  degree <- as.integer(table(c(node1, node2)))
  degree <- degree[degree > 0L]

  if (length(degree) == 0L) {
    return(data.frame(R2 = NA_real_))
  }

  degree_table <- table(degree)
  R2 <- if (length(unique(degree)) >= 3L) {
    fit <- lm(
      log10(as.numeric(degree_table) / sum(degree_table)) ~
        log10(as.numeric(names(degree_table)))
    )
    summary(fit)$r.squared
  } else {
    NA_real_
  }

  data.frame(R2 = R2)
}

string_links <- read.csv(
  "/proj/c.zihao/work1/01data/06string/links.csv",
  stringsAsFactors = FALSE
)

string_metrics <- bind_rows(lapply(c(0.05, 0.10, 0.15, 0.20), function(pct) {
  n_keep <- max(1L, round(pct * nrow(string_links)))
  selected <- order(string_links$score, decreasing = TRUE)[seq_len(n_keep)]
  metrics <- network_metrics(
    string_links$protein1[selected],
    string_links$protein2[selected]
  )
  data.frame(
    Method = "STRING",
    top_pct = pct,
    R2 = metrics$R2
  )
}))

network_metrics_table <- bind_rows(method_metrics, string_metrics)

method_order <- c(
  "STRING", "MOSSN",
  "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress",
  "noCorr", "noRWR", "noSeed",
  "RandomBackbone", "PermutedControl"
)
rule_order <- c("Method-specific (KS-optimal)", "STRING reference")

auto <- read.csv(file.path(result_dir, "gof_auto_summary.csv"), stringsAsFactors = FALSE) %>%
  transmute(
    Method        = method,
    top_pct,
    xmin_rule     = rule_order[1],
    x_min         = mean_xmin,
    KS            = mean_ks,
    Bootstrap_p   = mean_bootstrap_p,
    Frac_p_ge_0.1 = proportion_p_ge_0_1
  ) %>%
  left_join(network_metrics_table, by = c("Method", "top_pct"))

strg <- read.csv(file.path(result_dir, "gof_string_summary.csv"), stringsAsFactors = FALSE) %>%
  transmute(
    Method        = method,
    top_pct,
    xmin_rule     = rule_order[2],
    x_min         = string_xmin,
    KS            = mean_ks,
    Bootstrap_p   = mean_bootstrap_p,
    Frac_p_ge_0.1 = proportion_p_ge_0_1
  ) %>%
  left_join(network_metrics_table, by = c("Method", "top_pct"))

out <- bind_rows(auto, strg) %>%
  filter(Method %in% method_order) %>%
  mutate(
    Method    = factor(Method, levels = method_order),
    xmin_rule = factor(xmin_rule, levels = rule_order)
  ) %>%
  arrange(xmin_rule, Method, top_pct) %>%
  transmute(
    Method        = as.character(Method),
    Threshold     = paste0(top_pct * 100, "%"),
    xmin_rule     = as.character(xmin_rule),
    x_min         = sprintf("%.5f", x_min),
    R2            = sprintf("%.5f", R2),
    KS            = sprintf("%.5f", KS),
    Bootstrap_p   = sprintf("%.5f", Bootstrap_p),
    Frac_p_ge_0.1 = sprintf("%.5f", Frac_p_ge_0.1)
  )

names(out)[names(out) == "R2"] <- "R²"

out_file <- file.path(result_dir, "SupplementaryTable6_ScaleFreeGOF.csv")
write.csv(out, out_file, row.names = FALSE, quote = FALSE)
cat("Saved ->", out_file, "(", nrow(out), "rows )\n")
