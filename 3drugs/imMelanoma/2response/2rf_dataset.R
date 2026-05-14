# =================================================================
rm(list = ls())
library(randomForest)

# =================================================================

variableNumber <- 0.2
sampleNumber   <- 500
n_trees        <- 500
matrix_SeleNum <- c(1, 2, 3, 4, 5, 6, 7, 8)

dieasename <- "imMelanoma"
base_path  <- "/proj/c.zihao/work1/3drugs/"
surDir     <- paste0(base_path, dieasename, "/2response/")
netDir     <- paste0(base_path, dieasename, "/2response/net/")

setwd(surDir)
survdata <- read.csv("survdata.csv", row.names = 1)

print(head(survdata))
print(table(survdata$dataset))
summary(survdata$Response)

# =================================================================

lm_screen_topN <- function(train_df, response_col, feature_cols, topN) {
  pvals <- sapply(feature_cols, function(g) {
    fml <- as.formula(paste0(response_col, " ~ `", g, "`"))
    tryCatch({
      fit <- lm(fml, data = train_df)
      s   <- summary(fit)$coefficients
      if (nrow(s) < 2) return(1.0)
      pval <- s[2, "Pr(>|t|)"]
      if (!is.finite(pval)) 1.0 else pval
    }, error = function(e) 1.0)
  })
  names(pvals) <- feature_cols
  ord <- order(pvals, na.last = TRUE)
  names(pvals)[ord][seq_len(min(topN, length(pvals)))]
}

safe_cor <- function(pred, true, method = "spearman") {
  valid <- which(is.finite(pred) & is.finite(true))
  if (length(valid) < 3) return(NA_real_)
  if (var(true[valid]) == 0) return(NA_real_)
  if (var(pred[valid]) == 0) return(0)
  suppressWarnings(cor(pred[valid], true[valid], method = method))
}

# =================================================================

all_results <- data.frame(
  File        = character(),
  Dataset     = character(),
  Spearman    = numeric(),
  Pearson     = numeric(),
  N_test      = integer(),
  stringsAsFactors = FALSE
)

spearman_results <- numeric(length(matrix_SeleNum))

for (k in seq_along(matrix_SeleNum)) {
  i <- matrix_SeleNum[k]
  message(paste0("Processing PPPbi", i, ".csv"))

  net_path <- paste0(netDir, "PPPbi", i, ".csv")
  net <- read.csv(net_path, header = TRUE, row.names = 1)
  print(dim(net))

  net <- as.data.frame(t(net))
  net$Sample <- rownames(net)

  data <- merge(survdata, net, by = "Sample")
  rownames(data) <- data$Sample
  data$Sample <- NULL

  colnames(data) <- gsub("-", "_", colnames(data))

  data$Response <- as.numeric(data$Response)

  datasetnames <- unique(data$dataset)
  spearman_list <- c()

  for (ds in datasetnames) {

    print(ds)
    train_data <- data[data$dataset != ds, ]
    test_data  <- data[data$dataset == ds, ]

    train_data$dataset <- NULL
    test_data$dataset  <- NULL

    print(dim(train_data))
    if (nrow(train_data) > sampleNumber) {
      set.seed(1)
      train_data_subset <- train_data[sample(nrow(train_data), sampleNumber), ]
    } else {
      train_data_subset <- train_data
    }
    print(dim(train_data_subset))

    # Step 1: univariate linear regression screening
    feature_cols  <- setdiff(colnames(train_data_subset), "Response")
    top_feature_n <- max(1, ceiling(length(feature_cols) * variableNumber))

    top_features <- lm_screen_topN(
      train_df     = train_data_subset,
      response_col = "Response",
      feature_cols = feature_cols,
      topN         = top_feature_n
    )

    # Step 2: reduce to top features
    selected_features <- c("Response", top_features)
    train_data_sel <- train_data[, selected_features, drop = FALSE]
    test_data_sel  <- test_data[,  selected_features, drop = FALSE]

    # Step 3: random forest regression
    pred_score <- rep(mean(train_data_sel$Response, na.rm = TRUE),
                      nrow(test_data_sel))

    rf_fit <- tryCatch(
      suppressWarnings(
        randomForest(
          Response ~ .,
          data       = train_data_sel,
          ntree      = n_trees,
          importance = FALSE
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(rf_fit)) {
      pred_score <- as.numeric(
        predict(rf_fit, newdata = test_data_sel)
      )
    }

    # Step 4: evaluate
    true_resp  <- test_data_sel$Response
    valid_idx  <- which(is.finite(true_resp) & is.finite(pred_score))
    n_test     <- length(valid_idx)

    spearman_val <- safe_cor(pred_score, true_resp, method = "spearman")
    pearson_val  <- safe_cor(pred_score, true_resp, method = "pearson")

    spearman_list <- c(spearman_list, spearman_val)

    fmt_s <- ifelse(is.na(spearman_val), "NA", sprintf("%.4f", spearman_val))
    fmt_p <- ifelse(is.na(pearson_val),  "NA", sprintf("%.4f", pearson_val))
    message(sprintf("  %s | Spearman: %s | Pearson: %s | N: %d",
                    ds, fmt_s, fmt_p, n_test))

    all_results <- rbind(
      all_results,
      data.frame(
        File     = paste0("PPPbi", i),
        Dataset  = ds,
        Spearman = spearman_val,
        Pearson  = pearson_val,
        N_test   = n_test,
        stringsAsFactors = FALSE
      )
    )
  }

  print("===========================================================")
  spearman_results[k] <- mean(spearman_list, na.rm = TRUE)
}

final_results <- data.frame(
  File         = paste0("PPPbi", matrix_SeleNum),
  Mean_Spearman = round(spearman_results, 4)
)

print(final_results)
print(surDir)
setwd(surDir)
write.csv(all_results, "results_dataset.csv", row.names = FALSE)
