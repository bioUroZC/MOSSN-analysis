rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
source(paste0(PROJ_ROOT, "/function/LIONESS.R"))

base_dir    <- paste0(PROJ_ROOT, "/1NT/6time")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE)

n_reps <- 3

fixed_edges   <- 10000L
fixed_samples <- 10L

grid <- rbind(
  data.frame(n_samples = c(10L, 20L, 50L, 100L), n_edges = fixed_edges),
  data.frame(n_samples = fixed_samples, n_edges = c(20000L, 50000L, 100000L))
)

results <- data.frame(
  method            = character(),
  n_samples         = integer(),
  n_edges           = integer(),
  rep               = integer(),
  wall_time_s       = numeric(),
  time_per_sample_s = numeric(),
  peak_memory_mb    = numeric(),
  stringsAsFactors  = FALSE
)

for (i in seq_len(nrow(grid))) {
  N <- grid$n_samples[i]
  K <- grid$n_edges[i]

  expr_file  <- file.path(data_dir, sprintf("expr_S%03d.csv", N))
  links_file <- file.path(data_dir, sprintf("links_E%06d.csv", K))
  out_csv    <- file.path(results_dir, sprintf("timing_LIONESS_S%03d_E%06d.csv", N, K))

  for (rep in seq_len(n_reps) - 1L) {
    cat(sprintf("[LIONESS] N=%d E=%d rep=%d ...\n", N, K, rep))
    t_start <- proc.time()[["elapsed"]]
    tryCatch({
      LIONcal(exprSetFile = expr_file, ppiFile = links_file)
      wall <- proc.time()[["elapsed"]] - t_start
      cat(sprintf("  wall=%.2fs  /sample=%.4fs\n", wall, wall / N))
      results <- rbind(results, data.frame(
        method = "LIONESS", n_samples = N, n_edges = K, rep = rep,
        wall_time_s = round(wall, 4),
        time_per_sample_s = round(wall / N, 6),
        peak_memory_mb = NA_real_,
        stringsAsFactors = FALSE))
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
      results <<- rbind(results, data.frame(
        method = "LIONESS", n_samples = N, n_edges = K, rep = rep,
        wall_time_s = NA_real_, time_per_sample_s = NA_real_,
        peak_memory_mb = NA_real_,
        stringsAsFactors = FALSE))
    })
  }

  write.csv(results[results$n_samples == N & results$n_edges == K, ],
            out_csv, row.names = FALSE)
  cat(sprintf("  Saved → %s\n", out_csv))
}
