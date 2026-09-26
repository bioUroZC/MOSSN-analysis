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

clean_degree <- function(degree) {
  as.integer(degree[is.finite(degree) & degree > 0L])
}

fit_powerlaw_auto <- function(degree) {
  degree <- clean_degree(degree)
  if (length(degree) < 20L || length(unique(degree)) < 3L) return(NULL)
  tryCatch({
    model <- displ$new(degree)
    fit <- estimate_xmin(model, distance = "ks")
    list(xmin = as.integer(fit$xmin), n_tail = as.integer(fit$ntail))
  }, error = function(e) NULL)
}

fit_powerlaw_fixed_xmin <- function(degree, xmin) {
  degree <- clean_degree(degree)
  if (!is.finite(xmin) || xmin < 1L || length(degree) < 20L) return(NULL)
  tail_degree <- degree[degree >= xmin]
  if (length(tail_degree) < 20L || length(unique(tail_degree)) < 3L) return(NULL)

  tryCatch({
    model <- displ$new(degree)
    model$setXmin(as.integer(xmin))
    pars <- estimate_pars(model)$pars
    model$setPars(pars)
    list(
      xmin = as.integer(xmin),
      alpha = as.numeric(pars),
      ks = as.numeric(get_distance_statistic(model, distance = "ks")),
      n_tail = as.integer(length(tail_degree)),
      degree = degree
    )
  }, error = function(e) NULL)
}

bootstrap_ks_p_fixed_xmin <- function(fit, n_sims, seed) {
  set.seed(seed)
  degree <- fit$degree
  body <- degree[degree < fit$xmin]
  p_tail <- fit$n_tail / length(degree)
  simulated_ks <- rep(NA_real_, n_sims)

  for (b in seq_len(n_sims)) {
    simulated_ks[b] <- tryCatch({
      n_tail <- rbinom(1L, length(degree), p_tail)
      n_body <- length(degree) - n_tail
      simulated <- c(
        if (n_body > 0L && length(body) > 0L) sample(body, n_body, replace = TRUE) else integer(),
        if (n_tail > 0L) rpldis(n_tail, xmin = fit$xmin, alpha = fit$alpha) else integer()
      )
      simulated_fit <- fit_powerlaw_fixed_xmin(simulated, fit$xmin)
      if (is.null(simulated_fit)) NA_real_ else simulated_fit$ks
    }, error = function(e) NA_real_)
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
  dat <- read_csv(path, show_col_types = FALSE, col_select = all_of(c("Interaction", samples)))
  if (method == "RandomBackbone") dat <- dat[rowSums(dat[samples] != 0) > 0L, ]
  nodes <- strsplit(dat$Interaction, "_", fixed = TRUE)
  list(
    node1 = vapply(nodes, `[`, character(1), 1L),
    node2 = vapply(nodes, function(x) paste(x[-1L], collapse = "_"), character(1)),
    weights = dat[samples],
    samples = samples
  )
}

backbone <- read_csv(backbone_csv, show_col_types = FALSE) %>%
  select(score, protein1, protein2) %>%
  distinct()

backbone_degree <- function(top_pct) {
  selected <- edge_indices(backbone$score, top_pct, drop_zero = FALSE)
  as.integer(table(c(backbone$protein1[selected], backbone$protein2[selected])))
}

string_reference <- lapply(seq_along(top_pcts), function(pct_index) {
  degree <- backbone_degree(top_pcts[[pct_index]])
  auto <- fit_powerlaw_auto(degree)
  if (is.null(auto)) stop("Could not estimate a backbone xmin at top_pct ", top_pcts[[pct_index]])
  fit <- fit_powerlaw_fixed_xmin(degree, auto$xmin)
  list(xmin = auto$xmin, degree = degree, fit = fit)
})
string_xmins <- vapply(string_reference, `[[`, integer(1), "xmin")
cat("STRING backbone xmin per top_pct:",
    paste(sprintf("%.0f%%=%d", top_pcts * 100, string_xmins), collapse = ", "), "\n")

method_data <- setNames(lapply(methods, read_method_matrix), methods)
available_methods <- names(Filter(Negate(is.null), method_data))
if (length(available_methods) < 2L) stop("At least two method matrices are required.")
method_data <- method_data[available_methods]

common_samples <- Reduce(intersect, lapply(method_data, `[[`, "samples"))
common_samples <- intersect(luad_samples, common_samples)
if (max_samples > 0L && length(common_samples) > max_samples) {
  set.seed(20260915)
  common_samples <- sort(sample(common_samples, max_samples))
}
if (length(common_samples) == 0L) stop("No LUAD samples are shared by all available methods.")

degree_for_network <- function(dat, method, sample, top_pct) {
  selected <- edge_indices(
    as.numeric(dat$weights[[sample]]),
    top_pct = top_pct,
    drop_zero = !(method %in% signed_methods)
  )
  if (length(selected) == 0L) return(list(degree = integer(), n_edges = 0L, n_nodes = 0L))
  degree <- as.integer(table(c(dat$node1[selected], dat$node2[selected])))
  list(degree = degree, n_edges = length(selected), n_nodes = length(degree))
}

evaluate_at_string_xmin <- function(task_index) {
  sample_index <- ((task_index - 1L) %% length(common_samples)) + 1L
  pct_index <- ((task_index - 1L) %/% length(common_samples)) + 1L
  sample <- common_samples[[sample_index]]
  top_pct <- top_pcts[[pct_index]]
  string_xmin <- string_xmins[[pct_index]]

  bind_rows(lapply(seq_along(available_methods), function(method_index) {
    method <- available_methods[[method_index]]
    network <- degree_for_network(method_data[[method]], method, sample, top_pct)
    auto <- fit_powerlaw_auto(network$degree)
    fit <- fit_powerlaw_fixed_xmin(network$degree, string_xmin)
    seed <- pct_index * 10000000L + sample_index * 10000L + method_index

    data.frame(
      method = method,
      sample = sample,
      top_pct = top_pct,
      string_xmin = string_xmin,
      individual_xmin = if (is.null(auto)) NA_integer_ else auto$xmin,
      alpha = if (is.null(fit)) NA_real_ else fit$alpha,
      ks = if (is.null(fit)) NA_real_ else fit$ks,
      bootstrap_p = if (is.null(fit)) NA_real_ else bootstrap_ks_p_fixed_xmin(fit, n_boot, seed),
      n_tail = if (is.null(fit)) NA_integer_ else fit$n_tail,
      tail_fraction = if (is.null(fit) || network$n_nodes == 0L) NA_real_ else fit$n_tail / network$n_nodes,
      n_nodes = network$n_nodes,
      n_edges = network$n_edges,
      stringsAsFactors = FALSE
    )
  }))
}

task_count <- length(common_samples) * length(top_pcts)
results <- bind_rows(mclapply(seq_len(task_count), evaluate_at_string_xmin, mc.cores = n_cores))

backbone_results <- bind_rows(lapply(seq_along(top_pcts), function(pct_index) {
  ref <- string_reference[[pct_index]]
  fit <- ref$fit
  data.frame(
    method = "STRING",
    sample = "reference",
    top_pct = top_pcts[[pct_index]],
    string_xmin = ref$xmin,
    individual_xmin = ref$xmin,
    alpha = if (is.null(fit)) NA_real_ else fit$alpha,
    ks = if (is.null(fit)) NA_real_ else fit$ks,
    bootstrap_p = if (is.null(fit)) NA_real_ else bootstrap_ks_p_fixed_xmin(fit, n_boot, 9000000L + pct_index),
    n_tail = if (is.null(fit)) NA_integer_ else fit$n_tail,
    tail_fraction = if (is.null(fit)) NA_real_ else fit$n_tail / length(ref$degree),
    n_nodes = length(ref$degree),
    n_edges = NA_integer_,
    stringsAsFactors = FALSE
  )
}))

results <- bind_rows(backbone_results, results)

summary <- results %>%
  group_by(method, top_pct) %>%
  summarise(
    n_networks = n(),
    n_valid_fit = sum(is.finite(alpha)),
    string_xmin = first(string_xmin),
    mean_individual_xmin = mean(individual_xmin, na.rm = TRUE),
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_ks = mean(ks, na.rm = TRUE),
    mean_bootstrap_p = mean(bootstrap_p, na.rm = TRUE),
    proportion_p_ge_0_1 = mean(bootstrap_p >= 0.1, na.rm = TRUE),
    mean_n_tail = mean(n_tail, na.rm = TRUE),
    mean_tail_fraction = mean(tail_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(mean_alpha, mean_ks, mean_bootstrap_p,
                  proportion_p_ge_0_1, mean_tail_fraction),
                ~ round(.x, 5)))

write_csv(results, file.path(out_dir, "gof_string.csv"))
write_csv(summary, file.path(out_dir, "gof_string_summary.csv"))
cat("Saved scalefree_gof_string_xmin.csv and scalefree_gof_string_xmin_summary.csv\n")
