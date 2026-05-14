# =================================================================
rm(list = ls())
library(dplyr)
library(igraph)
library(tidyr)
library(data.table)
library(survival)
library(randomForestSRC)
library(survcomp)
library(timeROC)
library(caret)
library(pROC)

# =================================================================

dieasename <- "OV"

variableNumber <- 0.2
sampleNumber <- 500
ntreeNumber <- 500
matrix_SeleNum <- c(1, 2, 3, 4, 5, 6, 7, 8)
base_path <- "/proj/c.zihao/work1/2survival/"
surDir <- paste0(base_path, dieasename, "/3sur/")
netDir <- paste0(base_path, dieasename, "/3sur/net/")

setwd(surDir)
survdata <- read.csv('survdata.csv', row.names = 1)

head(survdata)
table(survdata$dataset)

survdata <- subset(survdata, survdata$OS_Time > 0)

cox_screen_topN <- function(train_df, time_col, event_col, feature_cols, topN = variableNumber) {
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
  finite_scores <- risk_scores[is.finite(risk_scores)]
  if (length(finite_scores) == 0 || length(unique(finite_scores)) <= 1) {
    return(0.5)
  }

  cidx <- tryCatch(
    concordance.index(
      x = risk_scores,
      surv.time = surv_time,
      surv.event = surv_event
    )$c.index,
    error = function(e) NA_real_
  )

  if (is.na(cidx)) 0.5 else cidx
}

safe_mean_auc <- function(risk_scores, surv_time, surv_event) {
  finite_scores <- risk_scores[is.finite(risk_scores)]
  if (length(finite_scores) == 0 || length(unique(finite_scores)) <= 1) {
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

  auc_val <- mean(td_auc$AUC, na.rm = TRUE)
  if (is.na(auc_val)) 0.5 else auc_val
}


# =================================================================

# =================================================================

# Prepare result storage
cindex_results <- numeric(length(matrix_SeleNum))
auc_results    <- numeric(length(matrix_SeleNum))

# Prepare result storage for each file and dataset
all_results <- data.frame(
  File = character(),
  Dataset = character(),
  C_index = numeric(),
  Mean_tAUC = numeric(),
  stringsAsFactors = FALSE
)


# Loop through PPPbi1.csv to PPPbi8.csv
for (k in seq_along(matrix_SeleNum)) {
   i <- matrix_SeleNum[k]
  message(paste("Processing PPPbi", i, ".csv", sep = ""))
  
  net_path <- paste0(netDir, "PPPbi", i, ".csv")
  net <- read.csv(net_path, header = TRUE, row.names = 1)
  print(dim(net))

  #if (nrow(net) > 5000) {
  #  cv <- apply(net, 1, function(x) sd(x) / mean(x))
  #  top_idx <- order(cv, decreasing = TRUE)[1:5000]
  #  net <- net[top_idx, ]
  #}
  
   
  net <- as.data.frame(t(net))
  net$Sample <- rownames(net)
  
  # Merge with survival data
  data <- merge(survdata, net, by = "Sample")
  rownames(data) <- data$Sample
  data$Sample <- NULL
  
  table(data$dataset)
  datasetnames <- unique(data$dataset)
  
  # Prepare per-dataset result storage
  cidx_list <- c()
  auc_list <- c()
  
  for (ds in datasetnames) {
    
    print(ds)
    train_data <- data[which(data$dataset != ds), ]
    test_data  <- data[which(data$dataset == ds), ]
    
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


    # Step 1: feature importance
    outcome_cols <- c("OS", "OS_Time")
    feature_cols <- setdiff(colnames(train_data_subset), outcome_cols)
    

    top_feature_n <- max(1, ceiling(length(feature_cols) * variableNumber))

    top_features <- cox_screen_topN(
      train_df = train_data_subset,     # IMPORTANT: use your subset object
      time_col = "OS_Time",
      event_col = "OS",
      feature_cols = feature_cols,
      topN = top_feature_n
    )
    
    selected_features <- c(outcome_cols, top_features)
    train_data_sel <- train_data[, selected_features, drop = FALSE]
    test_data_sel  <- test_data[,  selected_features, drop = FALSE]

    # Step 2: Reduce to top features
    selected_features <- c("OS", "OS_Time", top_features)
    train_data_sel <- train_data[, selected_features]
    test_data_sel <- test_data[, selected_features]
    
    # Step 3: Retrain RSF on selected features
    rsf_model <- rfsrc(Surv(OS_Time, OS) ~ ., data = train_data_sel, 
                       importance = FALSE, ntree = ntreeNumber)
    
    # Step 4: Predict
    rsf_pred <- predict(rsf_model, newdata = test_data_sel)
    risk_scores <- rsf_pred$predicted  
    
    sampled_test_scores_df <- data.frame(
      OS = test_data_sel$OS,
      Time = test_data_sel$OS_Time,
      risk_scores = risk_scores
    )
    
    # Step 5: C-index
    cidx <- safe_concordance_index(
      risk_scores = sampled_test_scores_df$risk_scores,
      surv_time = sampled_test_scores_df$Time,
      surv_event = sampled_test_scores_df$OS
    )
    cidx_list <- c(cidx_list, cidx)
    
    # Step 6: Time-dependent AUC
    auc_val <- safe_mean_auc(
      risk_scores = sampled_test_scores_df$risk_scores,
      surv_time = sampled_test_scores_df$Time,
      surv_event = sampled_test_scores_df$OS
    )
    auc_list <- c(auc_list, auc_val)
    
    # Save the result for this file and dataset!
    all_results <- rbind(
      all_results,
      data.frame(
        File = paste0("PPPbi", i),
        Dataset = ds,
        C_index = cidx,
        Mean_tAUC = auc_val,
        stringsAsFactors = FALSE
      )
    )
    
    message(sprintf("  %s | C-index: %.4f | Mean AUC: %.4f", ds, cidx, auc_val))
  }
  print('===========================================================')
  # Store mean across datasets for each PPPbi
  cindex_results[k] <- mean(cidx_list, na.rm = TRUE)
  auc_results[k]    <- mean(auc_list, na.rm = TRUE)
}

# Combine and display results
final_results <- data.frame(
  File = paste0("PPPbi", matrix_SeleNum),
  C_index = round(cindex_results, 4),
  Mean_tAUC = round(auc_results, 4)
)

print(final_results)
print(surDir)
setwd(surDir)
write.csv(all_results, 'rsf_dataset.csv')
