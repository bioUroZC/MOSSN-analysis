rm(list = ls())
library(dplyr)
library(igraph)
library(tidyr)
library(data.table)
library(randomForest)
library(caret)
library(pROC)
library(survival)
library(glmnet)
library(survcomp)
library(timeROC)

dieasename <- "GBM"
base_path <- "/proj/c.zihao/work1/03survival/"

saveDir <- paste0(base_path, dieasename, "/03sur/")

load_clinical_data <- function() {
  load_one <- function(name) {
    df <- read.csv(paste0(base_path, dieasename, "/", name, "/data/pd.csv"),
                   header = TRUE, row.names = 1)
    
    if ("OS" %in% colnames(df)) {
      df <- subset(df, select = c("Sample", "OS", "OS_Time"))
    } else {
      df <- subset(df, select = c("Sample", "EFS", "EFS_Time"))
      colnames(df)[2:3] <- c("OS", "OS_Time")
    }
    df$dataset <- name
    return(df)
  }
  datasets <- c("CGGA301", "CGGA325", "CGGA693", "GSE13041", "GSE16011",
                "GSE4412", "GSE72951", "TCGAGBM")
  
  do.call(rbind, lapply(datasets, load_one))
}

survdata <- load_clinical_data()

survdata$Sample <- make.names(survdata$Sample)
head(survdata)

getwd()

survdata <- subset(survdata, survdata$OS_Time > 0)
survdata <- survdata[!duplicated(survdata$Sample), ]
print('the number of samples in survival data')
print(dim(survdata))
setwd(saveDir)
write.csv(survdata, file = 'survdata.csv')

get_top_links <- function(expr_data, percent = NULL, top_k = NULL) {
    mean_expr <- rowMeans(expr_data, na.rm = TRUE)
    sd_expr <- apply(expr_data, 1, sd, na.rm = TRUE)
    cv_expr <- sd_expr / mean_expr
    
    df <- data.frame(
      link = rownames(expr_data),
      mean_expr = mean_expr,
      sd_expr = sd_expr,
      cv_expr = cv_expr,
      stringsAsFactors = FALSE
    )
    if (is.null(top_k)) {
      top_k <- ceiling(nrow(df) * percent)
    }
    top_k <- min(top_k, nrow(df))

    top_df <- df %>%
      dplyr::arrange(desc(cv_expr)) %>%
      dplyr::slice(1:top_k)
    
    return(top_df)
}

filter_shared_links <- function(expr_data, surv_df) {
    project_list <- unique(surv_df$dataset)

    keep_idx <- sapply(project_list, function(proj) {
      samples <- surv_df %>%
        dplyr::filter(dataset == proj) %>%
        dplyr::pull(Sample)

      sub_data <- expr_data[, colnames(expr_data) %in% samples, drop = FALSE]

      if (ncol(sub_data) == 0) {
        return(rep(FALSE, nrow(expr_data)))
      }

      rowSums(sub_data != 0, na.rm = TRUE) > 0
    })

    if (is.null(dim(keep_idx))) {
      keep_idx <- matrix(keep_idx, ncol = 1)
    }

    keep_links <- rowSums(keep_idx) == length(project_list)
    expr_data[keep_links, , drop = FALSE]
}

cox_screen_topN <- function(train_df, time_col, event_col, feature_cols, topN) {
  pvals <- sapply(feature_cols, function(g) {
    fml <- as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ `", g, "`"))
    tryCatch({
      fit <- coxph(fml, data = train_df, ties = "efron")
      s <- summary(fit)
      if (nrow(s$coefficients) == 0) return(1.0)
      pval <- s$coefficients[1, "Pr(>|z|)"]
      if (!is.finite(pval)) 1.0 else pval
    }, error = function(e) 1.0)
  })
  names(pvals) <- feature_cols
  ord <- order(pvals, na.last = TRUE)
  names(pvals)[ord][1:min(topN, length(pvals))]
}

safe_concordance_index <- function(risk_scores, surv_time, surv_event) {
  keep <- is.finite(risk_scores) & is.finite(surv_time) & is.finite(surv_event)
  risk_scores <- risk_scores[keep]
  surv_time <- surv_time[keep]
  surv_event <- surv_event[keep]

  if (length(unique(risk_scores)) <= 1) {
    return(0.5)
  }

  cidx <- tryCatch(
    concordance.index(
      x = risk_scores,
      surv.time = surv_time,
      surv.event = surv_event
    )$c.index,
    error = function(e) NA
  )

  if (is.na(cidx)) {
    return(0.5)
  }

  cidx
}

safe_mean_auc <- function(risk_scores, surv_time, surv_event) {
  keep <- is.finite(risk_scores) & is.finite(surv_time) & is.finite(surv_event)
  risk_scores <- risk_scores[keep]
  surv_time <- surv_time[keep]
  surv_event <- surv_event[keep]

  if (length(unique(risk_scores)) <= 1) {
    return(0.5)
  }

  timepoints <- quantile(surv_time, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  timepoints <- unique(as.numeric(timepoints))
  timepoints <- timepoints[is.finite(timepoints) & timepoints > 0]

  if (length(timepoints) == 0) {
    return(0.5)
  }

  td_auc <- tryCatch(
    timeROC(
      T = surv_time,
      delta = surv_event,
      marker = risk_scores,
      cause = 1,
      times = timepoints
    ),
    error = function(e) NULL
  )

  if (is.null(td_auc)) {
    return(0.5)
  }

  mean(td_auc$AUC, na.rm = TRUE)
}

n_stage1 <- 10000  
n_stage2 <- 1000   
n_stage3 <- 100    
sampleNumber <- 500

matrixfile <-  paste0("/proj/c.zihao/work1/03survival/", dieasename, "/02matrix")

matrix_files <- c("MOSSN.csv", "noRWR.csv",
                  "PPIXpress.csv", "Patkar.csv",
                  "RawExpr.csv", "NodeRWR.csv",
                  "noCorr.csv", "RandomBackbone.csv")

matrix_SeleNum <- seq_along(matrix_files)

cindex_results <- numeric(length(matrix_SeleNum))
auc_results    <- numeric(length(matrix_SeleNum))

all_results <- data.frame(
  File = character(),
  Dataset = character(),
  C_index = numeric(),
  Mean_tAUC = numeric(),
  stringsAsFactors = FALSE
)

for (k in seq_along(matrix_SeleNum)) {
  i <- matrix_SeleNum[k]
  file <- matrix_files[i]
  method_name <- sub("\\.csv$", "", file)
  message(paste("Processing", file, sep = " "))
  
  setwd(matrixfile)
  data0 <- fread(file) |> as.data.frame()
  rownames(data0) <- data0$Interaction
  data0$Interaction <- NULL
  colnames(data0) <- make.names(colnames(data0))
  
  print('Matrix dimensions before filtering:')
  print(dim(data0))
  
  data0 <- data0[,which(colnames(data0) %in% survdata$Sample), drop = FALSE]
  missing_samples <- setdiff(survdata$Sample, colnames(data0))
  if (length(missing_samples) > 0) {
    zero_cols <- as.data.frame(matrix(0, nrow=nrow(data0), ncol=length(missing_samples),
                                      dimnames=list(rownames(data0), missing_samples)))
    data0 <- cbind(data0, zero_cols)
  }
  
  print(data0[1:4, 1:4])
  print('original data')
  print(dim(data0))
  
  data0 <- abs(data0)
  colnames(data0) <- gsub("^Cancer_", "", colnames(data0))
  original_n_links <- nrow(data0)
  
  datasetnames <- unique(survdata$dataset)
  
  cidx_list <- c()
  auc_list <- c()
  
  for (ds in datasetnames) {
    
    print(ds)
    train_surv <- survdata[which(survdata$dataset != ds), ]
    
    data <- filter_shared_links(data0, train_surv)
    print('data after removing dataset-specific zero links')
    print(dim(data))
    
    train_cols <- colnames(data) %in% train_surv$Sample
    target_top_k <- n_stage1
    print('target number of top links (stage 1, fixed CV filter)')
    print(target_top_k)

    top_links <- get_top_links(data[, train_cols, drop = FALSE], top_k = target_top_k)
    selected_genes <- top_links$link
    
    print('top selected genes')
    print(head(selected_genes))
    print('number of selected genes:')
    print(length(selected_genes))
    
    data_filtered  <- data[which(rownames(data) %in% selected_genes), , drop = FALSE]
    print('filtered data')
    print(dim(data_filtered))
    
    if (nrow(data_filtered) < target_top_k) {
      set.seed(123)
      
      n_extra <- target_top_k - nrow(data_filtered)
      
      filtered_indices <- which(rownames(data0) %in% rownames(data_filtered))
      
      all_indices <- seq_len(original_n_links)
      remaining_indices <- setdiff(all_indices, filtered_indices)
      
      if (length(remaining_indices) == 0) {
        sampled_extra <- sample(all_indices, n_extra, replace = TRUE)
      } else if (length(remaining_indices) < n_extra) {
        sampled_extra <- sample(remaining_indices, n_extra, replace = TRUE)
      } else {
        sampled_extra <- sample(remaining_indices, n_extra, replace = FALSE)
      }
      
      data_filtered <- rbind(
        data_filtered,
        data0[sampled_extra, , drop = FALSE]
      )
    }
    
    PPP_binary <- data_filtered
    print('final number')
    print(dim(PPP_binary))
    
    net <- as.data.frame(t(PPP_binary))
    net$Sample <- rownames(net)
    
    data_ml <- merge(survdata, net, by = "Sample")
    rownames(data_ml) <- data_ml$Sample
    data_ml$Sample <- NULL
    
    table(data_ml$dataset)
    
    train_data <- data_ml[which(data_ml$dataset != ds), ]
    test_data  <- data_ml[which(data_ml$dataset == ds), ]
    
    train_data$dataset <- NULL
    test_data$dataset <- NULL
    
    print(dim(train_data))
    if (nrow(train_data) > sampleNumber) {
      set.seed(1)   
      train_data_subset <- train_data[sample(nrow(train_data), sampleNumber), ]
      } else {
        train_data_subset <- train_data
        }
    
    print(dim(train_data_subset))
    
    outcome_cols <- c("OS", "OS_Time")
    feature_cols <- setdiff(colnames(train_data_subset), outcome_cols)

    x_stage2 <- as.matrix(train_data_subset[, feature_cols, drop = FALSE])
    cor_vec <- suppressWarnings(
      cor(x_stage2, train_data_subset$OS_Time, use = "pairwise.complete.obs")[, 1]
    )
    cor_vec[!is.finite(cor_vec)] <- 0
    stage2_features <- names(sort(abs(cor_vec), decreasing = TRUE))[
      seq_len(min(n_stage2, length(cor_vec)))
    ]
    print('features remaining after stage 2 (correlation filter):')
    print(length(stage2_features))

    top_feature_n <- min(n_stage3, length(stage2_features))
    top_features <- cox_screen_topN(
      train_df = train_data_subset,
      time_col = "OS_Time",
      event_col = "OS",
      feature_cols = stage2_features,
      topN = top_feature_n
    )

    selected_features <- c(outcome_cols, top_features)
    train_data_sel <- train_data[, selected_features, drop = FALSE]
    test_data_sel  <- test_data[,  selected_features, drop = FALSE]
    
    x_train <- as.matrix(train_data_sel[, top_features, drop = FALSE])
    x_test <- as.matrix(test_data_sel[, top_features, drop = FALSE])
    y_train <- Surv(train_data_sel$OS_Time, train_data_sel$OS)
    
    train_var <- apply(x_train, 2, var, na.rm = TRUE)
    keep_features <- names(train_var)[is.finite(train_var) & train_var > 0]
    
    if (length(keep_features) > 0) {
      x_train <- x_train[, keep_features, drop = FALSE]
      x_test <- x_test[, keep_features, drop = FALSE]
    } else {
      x_train <- NULL
      x_test <- NULL
    }
    
    risk_scores <- rep(0, nrow(test_data_sel))
    
    if (!is.null(x_train) && ncol(x_train) > 0) {
      nfolds_use <- max(3, min(5, nrow(x_train)))
      
      glmnet_fit <- tryCatch({
        set.seed(42)
        foldid <- sample(rep(seq_len(nfolds_use), length.out = nrow(x_train)))
        cvfit <- cv.glmnet(
          x = x_train,
          y = y_train,
          family = "cox",
          alpha = 1,
          foldid = foldid
        )
        glmnet(
          x = x_train,
          y = y_train,
          family = "cox",
          alpha = 1,
          lambda = cvfit$lambda.min
        )
      }, error = function(e) NULL)
      
      if (!is.null(glmnet_fit)) {
        risk_scores <- as.numeric(
          predict(glmnet_fit, newx = x_test, type = "link")
        )
      } else if (length(keep_features) > 0) {
        fallback_feature <- keep_features[1]
        fallback_formula <- as.formula(
          paste0("Surv(OS_Time, OS) ~ `", fallback_feature, "`")
        )
        
        fallback_fit <- tryCatch(
          coxph(fallback_formula, data = train_data_sel, ties = "efron"),
          error = function(e) NULL
        )
        
        if (!is.null(fallback_fit)) {
          fallback_pred <- tryCatch(
            predict(fallback_fit, newdata = test_data_sel, type = "lp"),
            error = function(e) NULL
          )
          
          if (!is.null(fallback_pred)) {
            risk_scores <- as.numeric(fallback_pred)
          }
        }
      }
    }
    
    sampled_test_scores_df <- data.frame(
      OS = test_data_sel$OS,
      Time = test_data_sel$OS_Time,
      risk_scores = risk_scores
    )
    
    cidx <- safe_concordance_index(
      risk_scores = sampled_test_scores_df$risk_scores,
      surv_time = sampled_test_scores_df$Time,
      surv_event = sampled_test_scores_df$OS
    )
    cidx_list <- c(cidx_list, cidx)
    
    auc_val <- safe_mean_auc(
      risk_scores = sampled_test_scores_df$risk_scores,
      surv_time = sampled_test_scores_df$Time,
      surv_event = sampled_test_scores_df$OS
    )
    auc_list <- c(auc_list, auc_val)
    
    all_results <- rbind(
      all_results,
      data.frame(
        File = method_name,
        Dataset = ds,
        C_index = cidx,
        Mean_tAUC = auc_val,
        stringsAsFactors = FALSE
      )
    )
    
    message(sprintf("  %s | C-index: %.4f | Mean AUC: %.4f", ds, cidx, auc_val))
  }
  print('===========================================================')
  
  cindex_results[k] <- mean(cidx_list, na.rm = TRUE)
  auc_results[k]    <- mean(auc_list, na.rm = TRUE)
}

final_results <- data.frame(
  File = sub("\\.csv$", "", matrix_files[matrix_SeleNum]),
  C_index = round(cindex_results, 4),
  Mean_tAUC = round(auc_results, 4)
)

print(final_results)
print(saveDir)
setwd(saveDir)
all_results$C_index <- round(all_results$C_index, 4)
all_results$Mean_tAUC <- round(all_results$Mean_tAUC, 4)
write.csv(all_results, 'ml_dataset.csv')
