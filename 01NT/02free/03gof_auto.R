rm(list = ls())
library(dplyr)
library(readr)
library(parallel)
library(poweRlaw)

matrix_dir <- "/proj/c.zihao/work1/01NT/01string"
backbone_csv <- "/proj/c.zihao/work1/01data/06string/links.csv"
metadata_csv <- "/proj/c.zihao/work1/01data/03TCGA/metadata.csv"
out_dir <- "/proj/c.zihao/work1/01NT/02free/result"

methods <- c(
  "MOSSN", "noSeed", "noCorr", "noRWR", "RandomBackbone", "PermutedControl",
  "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress"
)
signed_methods <- c("SSN", "SWEET", "LIONESS")
top_pcts <- c(0.05, 0.10, 0.15, 0.20)
n_boot <- as.integer(Sys.getenv("SCALEFREE_BOOTSTRAPS", "200"))
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
max_samples <- as.integer(Sys.getenv("SCALEFREE_MAX_SAMPLES", "0"))

if (is.na(n_boot) || n_boot < 1L) stop("SCALEFREE_BOOTSTRAPS must be at least 1.")
if (is.na(max_samples) || max_samples < 0L) stop("SCALEFREE_MAX_SAMPLES must be non-negative.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_csv, stringsAsFactors = FALSE)
luad_samples <- metadata$Sample[metadata$Type == "LUAD"]

fit_powerlaw <- function(degree) {
  degree <- as.integer(degree[is.finite(degree) & degree > 0L])
  if (length(degree) < 20L || length(unique(degree)) < 3L) return(NULL)
  tryCatch({
    model <- displ$new(degree)
    xmin_fit <- estimate_xmin(model, distance = "ks")
    model$setXmin(xmin_fit$xmin)
    model$setPars(xmin_fit$pars)
    list(
      xmin = as.integer(xmin_fit$xmin),
      alpha = as.numeric(xmin_fit$pars),
      ks = as.numeric(xmin_fit$gof),
      n_tail = as.integer(xmin_fit$ntail),
      degree = degree
    )
  }, error = function(e) NULL)
}

bootstrap_ks_p <- function(fit, n_sims, seed) {
  set.seed(seed)
  degree <- fit$degree
  body <- degree[degree < fit$xmin]
  p_tail <- fit$n_tail / length(degree)
  simulated_ks <- rep(NA_real_, n_sims)

  for (b in seq_len(n_sims)) {
    n_tail <- rbinom(1L, length(degree), p_tail)
    n_body <- length(degree) - n_tail
    simulated <- c(
      if (n_body > 0L && length(body) > 0L) sample(body, n_body, replace = TRUE) else integer(),
      if (n_tail > 0L) rpldis(n_tail, xmin = fit$xmin, alpha = fit$alpha) else integer()
    )
    simulated_fit <- fit_powerlaw(simulated)
    if (!is.null(simulated_fit)) simulated_ks[b] <- simulated_fit$ks
  }

  simulated_ks <- simulated_ks[is.finite(simulated_ks)]
  if (length(simulated_ks) < ceiling(0.8 * n_sims)) return(NA_real_)
  (sum(simulated_ks >= fit$ks) + 1) / (length(simulated_ks) + 1)
}

edge_indices <- function(weight, top_pct, drop_zero) {
  nonzero <- if (drop_zero) which(weight != 0) else seq_along(weight)
  if (length(nonzero) == 0L) return(integer())
  if (drop_zero && round(0.20 * length(nonzero)) < 500L) {
    padding <- setdiff(seq_along(weight), nonzero)
    return(c(nonzero, sample(padding, max(0L, 500L - length(nonzero)))))
  }
  n_keep <- max(1L, round(top_pct * length(nonzero)))
  nonzero[order(weight[nonzero], decreasing = TRUE)[seq_len(n_keep)]]
}

evaluate_degree <- function(node1, node2, weight, top_pct, drop_zero, seed) {
  selected <- edge_indices(weight, top_pct, drop_zero)
  if (length(selected) == 0L) {
    return(data.frame(xmin = NA_integer_, alpha = NA_real_, ks = NA_real_,
                      bootstrap_p = NA_real_, n_tail = NA_integer_,
                      n_nodes = NA_integer_, n_edges = 0L))
  }
  degree <- as.integer(table(c(node1[selected], node2[selected])))
  fit <- fit_powerlaw(degree)
  if (is.null(fit)) {
    return(data.frame(xmin = NA_integer_, alpha = NA_real_, ks = NA_real_,
                      bootstrap_p = NA_real_, n_tail = NA_integer_,
                      n_nodes = length(degree), n_edges = length(selected)))
  }
  data.frame(
    xmin = fit$xmin,
    alpha = fit$alpha,
    ks = fit$ks,
    bootstrap_p = bootstrap_ks_p(fit, n_boot, seed),
    n_tail = fit$n_tail,
    n_nodes = length(degree),
    n_edges = length(selected)
  )
}

matrix_file <- function(method) {
  ablation_file <- file.path(matrix_dir, "01ablation", "02merge", "merged_matrices", method, "merged_matrix.csv")
  benchmark_file <- file.path(matrix_dir, "02benchmark", "02merge", "merged_matrices", method, "merged_matrix.csv")
  if (file.exists(ablation_file)) return(ablation_file)
  if (file.exists(benchmark_file)) return(benchmark_file)
  NA_character_
}

read_method_matrix <- function(method) {
  path <- matrix_file(method)
  if (is.na(path)) return(NULL)
  header <- names(read_csv(path, show_col_types = FALSE, n_max = 0))
  samples <- intersect(luad_samples, setdiff(header, "Interaction"))
  if (length(samples) == 0L) return(NULL)
  if (max_samples > 0L && length(samples) > max_samples) {
    set.seed(20260915)
    samples <- sort(sample(samples, max_samples))
  }
  dat <- read_csv(path, show_col_types = FALSE, col_select = all_of(c("Interaction", samples)))
  if (method == "RandomBackbone") dat <- dat[rowSums(dat[samples] != 0) > 0L, ]
  nodes <- strsplit(dat$Interaction, "_", fixed = TRUE)
  node1 <- vapply(nodes, `[`, character(1), 1L)
  node2 <- vapply(nodes, function(x) paste(x[-1L], collapse = "_"), character(1))
  list(node1 = node1, node2 = node2, weights = dat[samples], samples = samples)
}

process_method <- function(method_index, method) {
  cat("Processing", method, "\n")
  dat <- read_method_matrix(method)
  if (is.null(dat)) return(NULL)
  drop_zero <- !(method %in% signed_methods)
  task_grid <- expand.grid(sample_index = seq_along(dat$samples), pct_index = seq_along(top_pcts))
  rows <- mclapply(seq_len(nrow(task_grid)), function(task_index) {
    sample_index <- task_grid$sample_index[task_index]
    pct_index <- task_grid$pct_index[task_index]
    result <- evaluate_degree(
      dat$node1,
      dat$node2,
      as.numeric(dat$weights[[sample_index]]),
      top_pcts[pct_index],
      drop_zero,
      seed = method_index * 1000000L + sample_index * 1000L + pct_index
    )
    cbind(
      data.frame(method = method, sample = dat$samples[[sample_index]], top_pct = top_pcts[pct_index]),
      result
    )
  }, mc.cores = n_cores)
  bind_rows(rows)
}

backbone <- read_csv(backbone_csv, show_col_types = FALSE) %>%
  select(score, protein1, protein2) %>%
  distinct()
backbone_results <- bind_rows(lapply(seq_along(top_pcts), function(pct_index) {
  pct <- top_pcts[pct_index]
  n_keep <- max(1L, round(pct * nrow(backbone)))
  selected <- order(backbone$score, decreasing = TRUE)[seq_len(n_keep)]
  result <- evaluate_degree(
    backbone$protein1,
    backbone$protein2,
    backbone$score,
    pct,
    FALSE,
    seed = 9000000L + pct_index
  )
  cbind(data.frame(method = "STRING", sample = "reference", top_pct = pct), result)
}))

method_results <- bind_rows(Filter(
  Negate(is.null),
  Map(process_method, seq_along(methods), methods)
))
results <- bind_rows(backbone_results, method_results)

summary <- results %>%
  group_by(method, top_pct) %>%
  summarise(
    n_networks = n(),
    mean_xmin = mean(xmin, na.rm = TRUE),
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_ks = mean(ks, na.rm = TRUE),
    mean_bootstrap_p = mean(bootstrap_p, na.rm = TRUE),
    proportion_p_ge_0_1 = mean(bootstrap_p >= 0.1, na.rm = TRUE),
    mean_n_tail = mean(n_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(mean_alpha, mean_ks, mean_bootstrap_p, proportion_p_ge_0_1),
                ~ round(.x, 5)))

write_csv(results, file.path(out_dir, "gof_auto.csv"))
write_csv(summary, file.path(out_dir, "gof_auto_summary.csv"))
cat("Saved scalefree_gof.csv and scalefree_gof_summary.csv\n")
