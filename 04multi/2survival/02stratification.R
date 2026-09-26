rm(list = ls())
library(dplyr)
library(data.table)
library(survival)

RESULTS_DIR <- "/proj/c.zihao/work1/04multi/2survival/merge"
FILES_DIR   <- "/proj/c.zihao/work1/04multi/1prepare/06files"
OUT_DIR     <- "/proj/c.zihao/work1/04multi/2survival/result"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CANCERS   <- c("BLCA", "LIHC", "LUAD", "SARC", "STAD")
             
TOP_FRACS <- c(0.05, 0.10, 0.15, 0.20)

METHODS <- c(
    "EXP", "MET", "CNV",
    "NoCross", "Restart",
    "Direct", "NoDyn", "MultiLayer"
)

safe_logrank_p <- function(time, event, group) {
    if (length(unique(group)) < 2) return(NA_real_)
    fit <- tryCatch(
        survdiff(Surv(time, event) ~ group),
        error = function(e) NULL
    )
    if (is.null(fit) || length(fit$n) < 2) return(NA_real_)
    pchisq(fit$chisq, df = length(fit$n) - 1, lower.tail = FALSE)
}

orient_groups_by_survival <- function(time, event, group) {
    df <- data.frame(OSTime = time, OS = event, Group = factor(group))
    fit <- tryCatch(
        survfit(Surv(OSTime, OS) ~ Group, data = df),
        error = function(e) NULL
    )

    levs <- levels(df$Group)
    if (length(levs) != 2 || is.null(fit)) {
        out <- factor(ifelse(as.character(df$Group) == levs[1], "Good", "Poor"),
                      levels = c("Good", "Poor"))
        return(list(group = out, stats = c(Good = NA_real_, Poor = NA_real_)))
    }

    tbl <- summary(fit)$table
    if (is.null(dim(tbl))) {
        tbl <- matrix(tbl, nrow = 1, dimnames = list(names(summary(fit)$table)[1], names(summary(fit)$table)))
    }

    score_for <- function(rowname) {
        row <- tbl[rowname, ]
        med <- suppressWarnings(as.numeric(row["median"]))
        if (is.finite(med)) return(med)
        rmean <- suppressWarnings(as.numeric(row["rmean"]))
        if (is.finite(rmean)) return(rmean)
        NA_real_
    }

    row1 <- paste0("Group=", levs[1])
    row2 <- paste0("Group=", levs[2])
    s1 <- score_for(row1)
    s2 <- score_for(row2)

    poor_level <- if (is.na(s1) || is.na(s2)) {
        if (mean(event[group == levs[1]], na.rm = TRUE) >= mean(event[group == levs[2]], na.rm = TRUE)) levs[1] else levs[2]
    } else if (s1 <= s2) {
        levs[1]
    } else {
        levs[2]
    }
    good_level <- setdiff(levs, poor_level)

    out <- factor(ifelse(as.character(df$Group) == poor_level, "Poor", "Good"),
                  levels = c("Good", "Poor"))
    stats <- c(
        Good = if (good_level == levs[1]) s1 else s2,
        Poor = if (poor_level == levs[1]) s1 else s2
    )
    list(group = out, stats = stats)
}

safe_cox_metrics <- function(time, event, group) {
    if (length(unique(group)) < 2) {
        return(list(hr = NA_real_, hr_low = NA_real_, hr_high = NA_real_, p = NA_real_))
    }

    df <- data.frame(
        OSTime = time,
        OS = event,
        Group = factor(group, levels = c("Good", "Poor"))
    )

    fit <- tryCatch(
        coxph(Surv(OSTime, OS) ~ Group, data = df, ties = "efron"),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(list(hr = NA_real_, hr_low = NA_real_, hr_high = NA_real_, p = NA_real_))
    }

    s <- summary(fit)
    coef_row <- s$coefficients[1, , drop = FALSE]
    conf_row <- s$conf.int[1, , drop = FALSE]

    list(
        hr = unname(conf_row[1, "exp(coef)"]),
        hr_low = unname(conf_row[1, "lower .95"]),
        hr_high = unname(conf_row[1, "upper .95"]),
        p = unname(coef_row[1, "Pr(>|z|)"])
    )
}

summarise_mode <- function(df, hr_col) {
    df |>
        group_by(Cancer, Method) |>
        summarise(
            Best_FeatureFrac = FeatureFrac[which.max(NegLog10_LogRank)],
            Best_NegLog10_LogRank = max(NegLog10_LogRank, na.rm = TRUE),
            Best_LogRank_p = LogRank_p[which.max(NegLog10_LogRank)],
            Best_HR_Poor_vs_Good = .data[[hr_col]][which.max(NegLog10_LogRank)],
            .groups = "drop"
        ) |>
        arrange(Cancer, desc(Best_NegLog10_LogRank))
}

surv_map <- setNames(lapply(CANCERS, function(cancer) {
    df <- read.csv(file.path(FILES_DIR, paste0(cancer, "_OS.csv")),
                   row.names = 1, stringsAsFactors = FALSE)
    subset(df, OSTime > 0 & Type == cancer)
}), CANCERS)

balanced_results <- data.frame()
balanced_assignments <- data.frame()

for (method in METHODS) {
    cat("\n[", method, "]\n", sep = "")

    mat_file <- file.path(RESULTS_DIR, method, "merged_matrix.csv")
    if (!file.exists(mat_file)) {
        cat("  merged_matrix.csv not found, skip\n")
        next
    }

    mat <- fread(mat_file) |> as.data.frame()
    rownames(mat) <- mat$Interaction
    mat$Interaction <- NULL
    mat <- abs(mat)

    for (cancer in CANCERS) {
        surv_df <- surv_map[[cancer]]
        keep_samples <- intersect(surv_df$Sample, colnames(mat))
        if (length(keep_samples) < 10) {
            cat("  ", cancer, ": too few samples, skip\n", sep = "")
            next
        }

        sub_mat <- mat[, keep_samples, drop = FALSE]
        sub_mat[is.na(sub_mat)] <- 0

        means <- rowMeans(sub_mat)
        sds   <- apply(sub_mat, 1, sd)
        cv    <- sds / (abs(means) + 1e-9)
        cv    <- cv[is.finite(cv)]
        if (length(cv) == 0) next

        sample_df <- surv_df[match(keep_samples, surv_df$Sample), c("Sample", "OS", "OSTime")]

        for (top_frac in TOP_FRACS) {
            frac_label <- paste0(as.integer(top_frac * 100), "%")
            top_n <- max(2, floor(length(cv) * top_frac))
            keep_features <- names(sort(cv, decreasing = TRUE))[seq_len(min(top_n, length(cv)))]

            X_raw <- t(as.matrix(sub_mat[keep_features, keep_samples, drop = FALSE]))
            nonconst <- apply(X_raw, 2, sd) > 0
            X_raw <- X_raw[, nonconst, drop = FALSE]
            if (ncol(X_raw) < 2) next

            X <- scale(X_raw)

            pca_res <- prcomp(X, center = FALSE, scale. = FALSE)
            pc1 <- pca_res$x[, 1]
            split_cut <- median(pc1, na.rm = TRUE)
            split_groups <- ifelse(pc1 > split_cut, "High", "Low")
            if (sum(split_groups == "High") == 0 || sum(split_groups == "Low") == 0) {
                ord <- order(pc1)
                half <- floor(length(pc1) / 2)
                split_groups <- rep("High", length(pc1))
                split_groups[ord[seq_len(half)]] <- "Low"
            }
            split_groups <- factor(split_groups, levels = c("Low", "High"))
            oriented_bal <- orient_groups_by_survival(sample_df$OSTime, sample_df$OS, split_groups)
            surv_groups_bal <- oriented_bal$group
            logrank_p_bal <- safe_logrank_p(sample_df$OSTime, sample_df$OS, surv_groups_bal)
            cox_bal <- safe_cox_metrics(sample_df$OSTime, sample_df$OS, surv_groups_bal)
            if (!is.na(cox_bal$hr) && cox_bal$hr < 1) {
                surv_groups_bal <- factor(ifelse(surv_groups_bal == "Poor", "Good", "Poor"), levels = c("Good", "Poor"))
                oriented_bal$stats <- c(Good = oriented_bal$stats["Poor"], Poor = oriented_bal$stats["Good"])
                names(oriented_bal$stats) <- c("Good", "Poor")
                cox_bal <- safe_cox_metrics(sample_df$OSTime, sample_df$OS, surv_groups_bal)
            }

            balanced_results <- rbind(balanced_results, data.frame(
                Cancer           = cancer,
                Method           = method,
                FeatureFrac      = frac_label,
                N_samples        = nrow(X),
                N_features       = ncol(X),
                Good_n           = sum(surv_groups_bal == "Good"),
                Poor_n           = sum(surv_groups_bal == "Poor"),
                LogRank_p        = signif(logrank_p_bal, 4),
                NegLog10_LogRank = round(-log10(logrank_p_bal + 1e-300), 4),
                Cox_HR_Poor_vs_Good = round(cox_bal$hr, 4),
                Cox_HR_low95     = round(cox_bal$hr_low, 4),
                Cox_HR_high95    = round(cox_bal$hr_high, 4),
                Cox_p            = signif(cox_bal$p, 4),
                Median_Good      = round(oriented_bal$stats["Good"], 4),
                Median_Poor      = round(oriented_bal$stats["Poor"], 4)
            ))

            balanced_assignments <- rbind(balanced_assignments, data.frame(
                Cancer        = cancer,
                Method        = method,
                FeatureFrac   = frac_label,
                Sample        = rownames(X),
                PC1           = as.numeric(pc1),
                SplitGroup    = as.character(split_groups),
                SurvivalGroup = as.character(surv_groups_bal),
                OS            = sample_df$OS,
                OSTime        = sample_df$OSTime
            ))

            cat(sprintf(
                "  %s %s | balanced HR=%.3f p=%.4g\n",
                cancer, frac_label, cox_bal$hr, logrank_p_bal
            ))
        }
    }

    rm(mat)
    gc()
}

balanced_summary <- summarise_mode(balanced_results, "Cox_HR_Poor_vs_Good")

write.csv(balanced_results, file.path(OUT_DIR, "results.csv"), row.names = FALSE)
write.csv(balanced_assignments, file.path(OUT_DIR, "assignments.csv"), row.names = FALSE)
write.csv(balanced_summary, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)

cat("\nSaved ->", OUT_DIR, "\n")
