rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
})

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
surv_base <- "/proj/c.zihao/work1/2survival"
out_root <- file.path(base_dir, "case_survival")
case_root <- out_root

args <- commandArgs(trailingOnly = TRUE)

edge_name <- "MMP9_SPP1"
matrix_file_name <- "MOSSN_noPrior_Matrix.csv"
time_col <- "OS_Time"
event_col <- "OS"
sample_col <- "Sample"
default_cancers <- c(
  "ACC", "BLCA", "BRCA", "CHOL",
  "CRC", "GBM", "KIRC", "LGG",
  "LIHC", "LUAD", "OV", "PAAD",
  "PRAD", "STAD"
)

edge_arg <- args[grepl("^--edge=", args)]
if (length(edge_arg) > 0) {
  edge_name <- sub("^--edge=", "", edge_arg[1])
}

cancer_arg <- args[grepl("^--cancers=", args)]

list_cancer_dirs <- function(root_dir) {
  dirs <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  keep <- vapply(dirs, function(x) {
    dir.exists(file.path(root_dir, x, "1Matrix")) &&
      dir.exists(file.path(root_dir, x, "3sur"))
  }, logical(1))
  dirs[keep]
}

if (length(cancer_arg) > 0) {
  cancer_dirs <- strsplit(sub("^--cancers=", "", cancer_arg[1]), ",", fixed = TRUE)[[1]]
  cancer_dirs <- cancer_dirs[cancer_dirs != ""]
} else {
  cancer_dirs <- intersect(default_cancers, list_cancer_dirs(surv_base))
}

edge_dir_name <- gsub("[^A-Za-z0-9_]+", "_", edge_name)
out_dir <- file.path(out_root, edge_dir_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

case_summary_path <- file.path(case_root, edge_dir_name, "edge_case_summary.csv")
case_summary <- if (file.exists(case_summary_path)) {
  read.csv(case_summary_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}

read_matrix_file <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("link" %in% colnames(df)) {
    rownames(df) <- df$link
    df$link <- NULL
  } else if (ncol(df) >= 2) {
    rownames(df) <- df[[2]]
    df <- df[, -(1:2), drop = FALSE]
  }
  as.matrix(df)
}

read_surv_file <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  bad_names <- which(is.na(colnames(df)) | colnames(df) == "")
  if (length(bad_names) > 0) {
    df <- df[, -bad_names, drop = FALSE]
  }
  df
}

safe_cox <- function(df, score_col, time_col, event_col) {
  keep <- is.finite(df[[score_col]]) & is.finite(df[[time_col]]) & is.finite(df[[event_col]])
  df_fit <- df[keep, , drop = FALSE]
  
  if (nrow(df_fit) < 20) {
    return(NULL)
  }
  if (length(unique(df_fit[[event_col]])) < 2) {
    return(NULL)
  }
  if (stats::sd(df_fit[[score_col]], na.rm = TRUE) < 1e-8) {
    return(NULL)
  }
  
  fml <- stats::as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", score_col))
  fit <- tryCatch(
    survival::coxph(fml, data = df_fit, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  
  ss <- summary(fit)
  if (is.null(ss$coefficients) || nrow(ss$coefficients) == 0) {
    return(NULL)
  }
  
  data.frame(
    beta = unname(ss$coefficients[1, "coef"]),
    hr = unname(ss$coefficients[1, "exp(coef)"]),
    z = unname(ss$coefficients[1, "z"]),
    p_value = unname(ss$coefficients[1, "Pr(>|z|)"]),
    conf_low = unname(ss$conf.int[1, "lower .95"]),
    conf_high = unname(ss$conf.int[1, "upper .95"]),
    n_samples_used = nrow(df_fit),
    n_events_used = sum(df_fit[[event_col]] == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

meta_from_betas <- function(beta, se) {
  keep <- is.finite(beta) & is.finite(se) & se > 0
  beta <- beta[keep]
  se <- se[keep]
  if (length(beta) < 2) {
    return(NULL)
  }
  
  w <- 1 / (se ^ 2)
  beta_meta <- sum(w * beta) / sum(w)
  se_meta <- sqrt(1 / sum(w))
  z_meta <- beta_meta / se_meta
  p_meta <- 2 * stats::pnorm(abs(z_meta), lower.tail = FALSE)
  
  data.frame(
    beta_meta = beta_meta,
    hr_meta = exp(beta_meta),
    se_meta = se_meta,
    z_meta = z_meta,
    p_meta = p_meta,
    conf_low = exp(beta_meta - 1.96 * se_meta),
    conf_high = exp(beta_meta + 1.96 * se_meta),
    n_studies = length(beta),
    stringsAsFactors = FALSE
  )
}

resolve_edge_name <- function(edge_label, row_ids) {
  if (edge_label %in% row_ids) {
    return(edge_label)
  }
  
  parts <- strsplit(edge_label, "_", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    rev_label <- paste(rev(parts), collapse = "_")
    if (rev_label %in% row_ids) {
      return(rev_label)
    }
  }
  
  NA_character_
}

score_edge_in_matrix <- function(mat, edge_label) {
  matched_edge <- resolve_edge_name(edge_label, rownames(mat))
  if (is.na(matched_edge)) {
    return(NULL)
  }
  
  score <- as.numeric(mat[matched_edge, ])
  data.frame(
    Sample = colnames(mat),
    edge_query = edge_label,
    edge_matched = matched_edge,
    edge_score = score,
    stringsAsFactors = FALSE
  )
}

all_score_rows <- list()
all_cox_rows <- list()
presence_rows <- list()
dataset_filter_rows <- list()

for (cancer_name in cancer_dirs) {
  matrix_path <- file.path(surv_base, cancer_name, "1Matrix", matrix_file_name)
  surv_path <- file.path(surv_base, cancer_name, "3sur", "survdata.csv")
  
  if (!file.exists(matrix_path) || !file.exists(surv_path)) {
    next
  }
  
  message("Processing ", cancer_name)
  
  mat <- read_matrix_file(matrix_path)
  surv_df <- read_surv_file(surv_path)
  
  matched_edge <- resolve_edge_name(edge_name, rownames(mat))
  presence_rows[[length(presence_rows) + 1]] <- data.frame(
    cancer = cancer_name,
    edge_query = edge_name,
    edge_matched = matched_edge,
    edge_found = !is.na(matched_edge),
    stringsAsFactors = FALSE
  )
  
  if (is.na(matched_edge)) {
    next
  }
  
  if (!(sample_col %in% colnames(surv_df))) {
    next
  }
  
  surv_df <- surv_df %>%
    dplyr::filter(
      is.finite(.data[[time_col]]),
      is.finite(.data[[event_col]]),
      .data[[time_col]] > 0
    )
  
  common_samples <- intersect(colnames(mat), surv_df[[sample_col]])
  if (length(common_samples) < 20) {
    next
  }
  
  mat <- mat[, common_samples, drop = FALSE]
  surv_df <- surv_df %>% dplyr::filter(.data[[sample_col]] %in% common_samples)
  
  score_df <- score_edge_in_matrix(mat, edge_name)
  if (is.null(score_df) || nrow(score_df) == 0) {
    next
  }
  
  score_df <- score_df %>%
    dplyr::filter(Sample %in% common_samples) %>%
    dplyr::mutate(cancer = cancer_name)
  
  all_score_rows[[cancer_name]] <- score_df
  
  dataset_list <- unique(surv_df$dataset)
  
  for (dataset_name in dataset_list) {
    surv_one <- surv_df %>% dplyr::filter(dataset == dataset_name)
    merge_df <- surv_one %>%
      dplyr::left_join(score_df, by = c("Sample" = "Sample"))
    
    score_keep <- is.finite(merge_df$edge_score)
    score_vals <- merge_df$edge_score[score_keep]
    dataset_all_zero <- length(score_vals) > 0 && all(score_vals == 0)
    
    dataset_filter_rows[[length(dataset_filter_rows) + 1]] <- data.frame(
      cancer = cancer_name,
      dataset = dataset_name,
      edge_query = edge_name,
      edge_matched = matched_edge,
      n_samples = nrow(merge_df),
      n_samples_with_score = length(score_vals),
      all_zero_edge_score = dataset_all_zero,
      stringsAsFactors = FALSE
    )
    
    if (dataset_all_zero) {
      message("Skipping ", cancer_name, " / ", dataset_name, ": edge score is all zero.")
      next
    }
    
    fit_df <- safe_cox(
      df = merge_df,
      score_col = "edge_score",
      time_col = time_col,
      event_col = event_col
    )
    
    if (is.null(fit_df)) {
      next
    }
    
    fit_df <- fit_df %>%
      dplyr::mutate(
        cancer = cancer_name,
        dataset = dataset_name,
        edge_query = edge_name,
        edge_matched = matched_edge,
        n_samples = nrow(merge_df),
        n_events = sum(merge_df[[event_col]] == 1, na.rm = TRUE),
        median_score = stats::median(merge_df$edge_score, na.rm = TRUE)
      )
    
    all_cox_rows[[length(all_cox_rows) + 1]] <- fit_df
  }
}

presence_tbl <- dplyr::bind_rows(presence_rows)
if (nrow(presence_tbl) > 0) {
  utils::write.csv(presence_tbl, file.path(out_dir, "edge_presence_by_cancer.csv"), row.names = FALSE)
}

dataset_filter_tbl <- dplyr::bind_rows(dataset_filter_rows)
if (nrow(dataset_filter_tbl) > 0) {
  utils::write.csv(dataset_filter_tbl, file.path(out_dir, "edge_dataset_filter_log.csv"), row.names = FALSE)
}

score_long <- dplyr::bind_rows(all_score_rows)
if (nrow(score_long) > 0) {
  utils::write.csv(score_long, file.path(out_dir, "edge_scores_survival_long.csv"), row.names = FALSE)
}

cox_tbl <- dplyr::bind_rows(all_cox_rows)
if (nrow(cox_tbl) == 0) {
  stop("No valid Cox results were generated for the requested edge.")
}

cox_tbl <- cox_tbl %>%
  dplyr::mutate(se = abs(beta / z)) %>%
  dplyr::arrange(p_value, cancer, dataset)

utils::write.csv(cox_tbl, file.path(out_dir, "edge_cox_per_dataset.csv"), row.names = FALSE)

meta_rows <- list()

meta_all <- meta_from_betas(cox_tbl$beta, cox_tbl$se)
if (!is.null(meta_all)) {
  meta_rows[[length(meta_rows) + 1]] <- cbind(
    edge_query = edge_name,
    scope = "all_datasets",
    meta_all,
    stringsAsFactors = FALSE
  )
}

for (cancer_name in unique(cox_tbl$cancer)) {
  df_ca <- cox_tbl %>% dplyr::filter(cancer == cancer_name)
  meta_ca <- meta_from_betas(df_ca$beta, df_ca$se)
  if (!is.null(meta_ca)) {
    meta_rows[[length(meta_rows) + 1]] <- cbind(
      edge_query = edge_name,
      scope = cancer_name,
      meta_ca,
      stringsAsFactors = FALSE
    )
  }
}

meta_tbl <- dplyr::bind_rows(meta_rows) %>%
  dplyr::mutate(
    edge_matched = edge_name,
    fdr = stats::p.adjust(p_meta, method = "BH")
  ) %>%
  dplyr::arrange(fdr, p_meta, dplyr::desc(abs(log(hr_meta))))

if (nrow(case_summary) > 0) {
  meta_tbl <- meta_tbl %>%
    dplyr::mutate(
      module_id = case_summary$module_id[1],
      module_direction = case_summary$direction[1],
      module_annotation = case_summary$concise_label[1],
      module_source = case_summary$best_source_hk[1]
    )
}

utils::write.csv(meta_tbl, file.path(out_dir, "edge_cox_meta_summary.csv"), row.names = FALSE)

plot_df <- meta_tbl %>%
  dplyr::filter(scope != "all_datasets") %>%
  dplyr::mutate(
    scope = factor(scope, levels = scope[order(hr_meta)]),
    direction = dplyr::if_else(hr_meta >= 1, "risk", "protective")
  )

if (nrow(plot_df) > 0) {
  p1 <- ggplot(plot_df, aes(x = hr_meta, y = scope, color = direction)) +
    geom_point(size = 2.8) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    scale_color_manual(values = c(risk = "#B2182B", protective = "#2166AC")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Meta hazard ratio by cancer",
      y = "Cancer type",
      color = "Direction",
      title = edge_name,
      subtitle = if (nrow(case_summary) > 0) {
        paste0(case_summary$module_id[1], " | ",
               case_summary$best_source_hk[1], " | ",
               case_summary$concise_label[1])
      } else {
        NULL
      }
    )
  
  pdf(file.path(out_dir, "edge_survival_meta_by_cancer.pdf"), height = 6, width = 8)
  print(p1)
  dev.off()
}

dataset_plot_df <- cox_tbl %>%
  dplyr::mutate(
    dataset_label = paste(cancer, dataset, sep = " | "),
    dataset_label = factor(dataset_label, levels = dataset_label[order(hr)]),
    direction = dplyr::if_else(hr >= 1, "risk", "protective")
  )

if (nrow(dataset_plot_df) > 0) {
  p2 <- ggplot(dataset_plot_df, aes(x = hr, y = dataset_label, color = direction)) +
    geom_point(size = 2.2) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.16) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    scale_color_manual(values = c(risk = "#B2182B", protective = "#2166AC")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Hazard ratio per dataset",
      y = "Cancer | Dataset",
      color = "Direction",
      title = paste0(edge_name, " prognosis across datasets")
    )
  
  pdf(file.path(out_dir, "edge_survival_per_dataset.pdf"), height = 10, width = 10)
  print(p2)
  dev.off()
}

message("Saved edge survival outputs to: ", out_dir)
