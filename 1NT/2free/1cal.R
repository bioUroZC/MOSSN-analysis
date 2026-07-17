rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(readr)
library(tidyr)
library(parallel)
library(poweRlaw)

# ── STRING-only config ────────────────────────────────────────────────────────

network <- "STRING"
matrix_dir   <- paste0(PROJ_ROOT, "/1NT/2string")
backbone_csv <- paste0(PROJ_ROOT, "/1NT/1data/string/links.csv")
# ablation methods live under matrix_dir/ablation/<method>, benchmark methods under matrix_dir/benchmark/<method>
methods      <- c(
  "MOSSN_uniform", "MOSSN_noSeed", "MOSSN_noCorr",
  "EdgeNoRWR", "DegreePreservedBackbone", "RandomBackbone", "PermutedControl",
  "SSN", "SWEET", "LIONESS", "Patkar", "PPIXpress", "Proteinarium"
)

n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
out_dir <- file.path(paste0(PROJ_ROOT, "/1NT/2free"), network)
dir.create(out_dir, showWarnings = FALSE)
top_pcts <- c(0.05, 0.10, 0.15, 0.20)

cat("Network:", network, "\n")
cat("Using", n_cores, "cores\n")

# ── LUAD sample list ──────────────────────────────────────────────────────────
luad_samples <- setdiff(
  names(read_csv(paste0(PROJ_ROOT, "/1NT/2string/benchmark/SSN/LUAD/Zscore.csv"),
                 show_col_types = FALSE, n_max = 0)),
  c("Gene1", "Gene2")
)
cat("LUAD samples:", length(luad_samples), "\n")

# ── backbone scale-free stats ─────────────────────────────────────────────────
backbone_stats <- local({
  df  <- read_csv(backbone_csv, show_col_types = FALSE) %>%
    select(protein1, protein2) %>% distinct()
  deg <- as.integer(table(c(df$protein1, df$protein2)))
  deg <- deg[deg > 0L]
  tbl <- table(deg)
  fit <- lm(log10(as.numeric(tbl) / sum(tbl)) ~ log10(as.numeric(names(tbl))))
  pl  <- tryCatch({
    o <- displ$new(deg); o$setXmin(estimate_xmin(o))
    o$setPars(estimate_pars(o)); o
  }, error = function(e) NULL)
  gamma_pl <- if (!is.null(pl)) pl$pars else NA_real_

  lr_stat <- NA_real_; lr_p <- NA_real_
  if (!is.null(pl)) tryCatch({
    ex <- disexp$new(deg); ex$setXmin(pl$getXmin())
    ex$setPars(estimate_pars(ex))
    comp <- compare_distributions(pl, ex)
    lr_stat <- comp$test_statistic; lr_p <- comp$p_two_sided
  }, error = function(e) NULL)

  list(R2 = round(summary(fit)$r.squared, 4),
       gamma_lm = round(-coef(fit)[[2L]], 4),
       gamma_pl = round(gamma_pl, 4),
       lr_stat  = round(lr_stat, 4),
       lr_p     = round(lr_p, 4),
       n_edges  = nrow(df),
       n_nodes  = length(deg))
})
cat(sprintf("%s backbone: %d edges, %d nodes, R2=%.4f, gamma_pl=%.4f, lr_stat=%.4f, lr_p=%.4f\n",
            network, backbone_stats$n_edges, backbone_stats$n_nodes,
            backbone_stats$R2, backbone_stats$gamma_pl,
            backbone_stats$lr_stat, backbone_stats$lr_p))

backbone_out <- as.data.frame(backbone_stats) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(backbone_out, file.path(out_dir, "backbone_stats.csv"))
cat("Backbone stats saved.\n")

# ── helpers ───────────────────────────────────────────────────────────────────

fit_scalefree <- function(n1, n2, w, top_pct) {
  empty <- list(R2 = NA_real_, gamma_lm = NA_real_, gamma_pl = NA_real_,
                lr_exp_stat = NA_real_, lr_exp_p = NA_real_,
                lr_ln_stat  = NA_real_, lr_ln_p  = NA_real_,
                H = NA_real_, degree_gini = NA_real_,
                top5pct_node_edge_frac = NA_real_)
  n <- length(w)
  if (n == 0L) return(empty)
  k   <- max(1L, round(top_pct * n))
  ord <- order(w, decreasing = TRUE)[seq_len(k)]
  deg_tbl <- table(c(n1[ord], n2[ord]))
  deg <- setNames(as.integer(deg_tbl), names(deg_tbl))
  deg <- deg[deg > 0L]
  if (length(unique(deg)) < 3L) return(empty)

  # degree heterogeneity H = <k^2> / <k>^2
  H <- mean(deg^2) / mean(deg)^2

  # degree Gini: 0 means uniform degree, 1 means highly concentrated hubs
  deg_sorted <- sort(deg)
  n_deg <- length(deg_sorted)
  degree_gini <- if (n_deg <= 1L || sum(deg_sorted) == 0) {
    0
  } else {
    (2 * sum(seq_len(n_deg) * deg_sorted)) / (n_deg * sum(deg_sorted)) -
      (n_deg + 1) / n_deg
  }

  # fraction of all edges touching the top 5% highest-degree nodes
  n_top_nodes <- max(1L, ceiling(0.05 * n_deg))
  top_nodes <- names(sort(deg, decreasing = TRUE))[seq_len(n_top_nodes)]
  edge_in_top <- (n1[ord] %in% top_nodes) | (n2[ord] %in% top_nodes)
  top5pct_node_edge_frac <- mean(edge_in_top)

  tbl    <- table(deg)
  fit_lm <- lm(log10(as.numeric(tbl) / sum(tbl)) ~ log10(as.numeric(names(tbl))))
  R2       <- summary(fit_lm)$r.squared
  gamma_lm <- -coef(fit_lm)[[2L]]

  # MLE power-law fit
  pl <- tryCatch({
    o <- displ$new(deg); o$setXmin(estimate_xmin(o))
    o$setPars(estimate_pars(o)); o
  }, error = function(e) NULL)
  gamma_pl <- if (!is.null(pl)) pl$pars else NA_real_

  lr_exp_stat <- NA_real_; lr_exp_p <- NA_real_
  lr_ln_stat  <- NA_real_; lr_ln_p  <- NA_real_

  if (!is.null(pl)) {
    xmin <- pl$getXmin()

    # LR test: power-law vs exponential
    tryCatch({
      ex <- disexp$new(deg); ex$setXmin(xmin)
      ex$setPars(estimate_pars(ex))
      comp        <- compare_distributions(pl, ex)
      lr_exp_stat <- comp$test_statistic
      lr_exp_p    <- comp$p_two_sided
    }, error = function(e) NULL)

    # LR test: power-law vs log-normal
    tryCatch({
      ln <- dislnorm$new(deg); ln$setXmin(xmin)
      ln$setPars(estimate_pars(ln))
      comp       <- compare_distributions(pl, ln)
      lr_ln_stat <- comp$test_statistic
      lr_ln_p    <- comp$p_two_sided
    }, error = function(e) NULL)
  }

  list(R2 = R2, gamma_lm = gamma_lm, gamma_pl = gamma_pl,
       lr_exp_stat = lr_exp_stat, lr_exp_p = lr_exp_p,
       lr_ln_stat  = lr_ln_stat,  lr_ln_p  = lr_ln_p,
       H = H, degree_gini = degree_gini,
       top5pct_node_edge_frac = top5pct_node_edge_frac)
}

read_matrix_long <- function(method) {
  f <- file.path(matrix_dir, "ablation", method, "merged_matrix.csv")
  if (!file.exists(f)) f <- file.path(matrix_dir, "benchmark", method, "merged_matrix.csv")
  if (!file.exists(f)) return(NULL)
  header  <- names(read_csv(f, show_col_types = FALSE, n_max = 0))
  samples <- intersect(setdiff(header, "Interaction"), luad_samples)
  if (length(samples) == 0L) return(NULL)
  dat <- read_csv(f, show_col_types = FALSE,
                  col_select = all_of(c("Interaction", samples)))
  if (method %in% c("DegreePreservedBackbone", "RandomBackbone")) {
    # These methods' merged matrix spans independently-generated topologies
    # per cancer; rows that are all-zero across LUAD samples belong to
    # another cancer's topology and aren't part of this network at all.
    dat <- dat[rowSums(dat[samples] != 0) > 0L, ]
  }
  dat %>%
    pivot_longer(all_of(samples), names_to = "sample", values_to = "weight") %>%
    separate("Interaction", into = c("node1", "node2"), sep = "_", extra = "merge") %>%
    select("sample", "node1", "node2", "weight")
}

process_method <- function(method_name) {
  cat(method_name, "...\n")
  dat <- read_matrix_long(method_name)
  if (is.null(dat) || nrow(dat) == 0L) return(NULL)
  dat_split <- split(dat, dat$sample, drop = TRUE)
  bind_rows(lapply(top_pcts, function(pct) {
    bind_rows(lapply(dat_split, function(grp) {
      m <- fit_scalefree(grp$node1, grp$node2, grp$weight, pct)
      data.frame(method = method_name, sample = grp$sample[1L], top_pct = pct,
                 R2 = m$R2, gamma_lm = m$gamma_lm, gamma_pl = m$gamma_pl,
                 lr_exp_stat = m$lr_exp_stat, lr_exp_p = m$lr_exp_p,
                 lr_ln_stat  = m$lr_ln_stat,  lr_ln_p  = m$lr_ln_p,
                 H = m$H, degree_gini = m$degree_gini,
                 top5pct_node_edge_frac = m$top5pct_node_edge_frac)
    }))
  }))
}

# ── compute ───────────────────────────────────────────────────────────────────

results <- bind_rows(Filter(
  Negate(is.null),
  mclapply(setNames(methods, methods), process_method, mc.cores = n_cores)
))

results <- results %>%
  mutate(across(where(is.double), ~ round(.x, 3)))

write_csv(results, file.path(out_dir, "scalefree_persample.csv"))
cat("Per-sample results saved.\n")
cat("Done.\n")
