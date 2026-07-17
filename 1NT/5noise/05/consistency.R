rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  match_idx <- grep("--file=", cmd_args)
  if (length(match_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[match_idx[1]]), mustWork = TRUE))
  }
  frame_files <- vapply(sys.frames(), function(env) {
    if (exists("ofile", envir = env, inherits = FALSE))
      get("ofile", envir = env, inherits = FALSE)
    else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0)
    return(normalizePath(frame_files[length(frame_files)], mustWork = TRUE))
  stop("Unable to determine script path.")
}


load_matrix <- function(path) {
  read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}


resolve_matrix_file <- function(root_dir, method) {
  candidates <- c(
    file.path(root_dir, method, "merged_matrix.csv"),
    file.path(root_dir, method, "merged.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}


load_best_matrix <- function(root_dir, method) {
  matrix_file <- resolve_matrix_file(root_dir, method)
  if (is.na(matrix_file)) return(NULL)
  load_matrix(matrix_file)
}


base_dir    <- dirname(get_script_path())
full_dir    <- paste0(PROJ_ROOT, "/1NT/2string")
sample_file <- file.path(base_dir, "data", "LUAD_noise_samples.csv")
methods     <- c("MOSSN_uniform", "SSN", "SWEET",
                 "LIONESS", "Patkar", "PPIXpress", "Proteinarium")
out_dir     <- file.path(base_dir, "consistency")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_samples <- read.csv(sample_file, stringsAsFactors = FALSE)$sample

consistency_rows <- list()
sample_rows      <- list()

for (method in methods) {
  sub_df  <- load_best_matrix(base_dir, method)
  full_df <- load_best_matrix(full_dir, method)

  if (is.null(sub_df) || is.null(full_df)) {
    cat(method, ": skipped (matrix not found)\n")
    next
  }

  common_edges   <- intersect(rownames(sub_df), rownames(full_df))
  common_samples <- selected_samples[
    selected_samples %in% colnames(sub_df) &
    selected_samples %in% colnames(full_df)
  ]

  if (length(common_edges) == 0 || length(common_samples) == 0) {
    cat(method, ": skipped (no common edges or samples)\n")
    next
  }

  sub_df  <- sub_df[common_edges, common_samples, drop = FALSE]
  full_df <- full_df[common_edges, common_samples, drop = FALSE]

  sample_spearman <- vapply(common_samples, function(s) {
    suppressWarnings(cor(as.numeric(full_df[[s]]), as.numeric(sub_df[[s]]),
                         method = "spearman", use = "pairwise.complete.obs"))
  }, numeric(1))

  sample_rows[[length(sample_rows) + 1]] <- data.frame(
    method = method, sample = common_samples,
    spearman = sample_spearman, stringsAsFactors = FALSE
  )

  consistency_rows[[length(consistency_rows) + 1]] <- data.frame(
    method = method,
    consistency = mean(sample_spearman, na.rm = TRUE),
    mean_spearman = mean(sample_spearman, na.rm = TRUE),
    n_samples = length(common_samples),
    stringsAsFactors = FALSE
  )

  cat(method, ": n=", length(common_samples),
      " spearman=", round(mean(sample_spearman, na.rm = TRUE), 4), "\n")
}

consistency_df <- if (length(consistency_rows) > 0) {
  out <- do.call(rbind, consistency_rows)
  out[order(-out$consistency), ]
} else {
  data.frame()
}

sample_df <- if (length(sample_rows) > 0) do.call(rbind, sample_rows) else data.frame()

write.csv(consistency_df, file.path(out_dir, "consistency_df.csv"), row.names = FALSE)
write.csv(sample_df, file.path(out_dir, "sample_level_consistency.csv"), row.names = FALSE)
