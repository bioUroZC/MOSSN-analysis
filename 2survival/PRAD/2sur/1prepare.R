# =================================================================
#                           
# =================================================================

rm(list = ls())
library(dplyr)
library(igraph)
library(tidyr)
library(data.table)
library(randomForest)
library(caret)
library(pROC)

# =================================================================

# =================================================================

dieasename <- "PRAD"
base_path <- "/proj/c.zihao/work1/2survival/"

saveDir <- paste0(base_path, dieasename, "/2sur/")
netDir <- paste0(base_path, dieasename, "/2sur/net/")
dir.create(netDir, recursive = TRUE, showWarnings = FALSE)

# =================================================================

# =================================================================

load_clinical_data <- function() {
  load_one <- function(name) {
    df <- read.csv(paste0(base_path, dieasename, "/", name, "/data/pd.csv"),
                   header = TRUE, row.names = 1)
    df <- subset(df, select = c("Sample", "OS", "OS_Time"))
    df$dataset <- name
    return(df)
  }
  datasets <- c("DKFZ2018", "GSE116918", "GSE21034", "GSE46602",
                "GSE54460", "GSE70768", "GSE70769", "TCGAPRAD")
  
  do.call(rbind, lapply(datasets, load_one))
}

survdata <- load_clinical_data()
head(survdata)

getwd()

survdata <- subset(survdata, survdata$OS_Time > 0)
print('the number of samples in survival data')
print(dim(survdata))
setwd(saveDir)
write.csv(survdata, file = 'survdata.csv')

# =================================================================

# =================================================================


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
  

matrixfile <-  paste0("/proj/c.zihao/work1/2survival/", dieasename, "/1Matrix")

matrix_files <- c("MOSSN_noPrior.csv", "MOSSN_uniform.csv",
                  "Patkar.csv", "PPIXpress.csv",
                  "Proteinarium.csv", "SSN.csv",
                  "SWEET.csv", "LIONESS.csv")

result_dims <- data.frame(
  matrix_file = character(),
  nrow = integer(),
  ncol = integer(),
  stringsAsFactors = FALSE
)


for (i in 1:8) {
  
  fileNum <- i
  
  file <- matrix_files[fileNum]
  cat("Processing file:", file, "\n")
  setwd(matrixfile)
  
  data <- fread(file) |> as.data.frame()
  rownames(data) <- data$Interaction
  data$Interaction <- NULL
  

  print('Matrix dimensions before filtering:')
  print(dim(data))
  
  
  data <- data[,which(colnames(data) %in% survdata$Sample)]
  missing_samples <- setdiff(survdata$Sample, colnames(data))
  if (length(missing_samples) > 0) {
    zero_cols <- as.data.frame(matrix(0, nrow=nrow(data), ncol=length(missing_samples),
                                      dimnames=list(rownames(data), missing_samples)))
    data <- cbind(data, zero_cols)
  }

  print(data[1:4, 1:4])
  print('original data')
  print(dim(data))
  
  data <- abs(data)
  
  colnames(data) <- gsub("^Cancer_", "", colnames(data))
  data_all <- data

  original_n_links <- nrow(data)

  data <- filter_shared_links(data, survdata)
  print('data after removing dataset-specific zero links')
  print(dim(data))
  
  percent <- 0.05
  print(percent)
  target_top_k <- ceiling(original_n_links * percent)
  print('target number of top links based on original edge count')
  print(target_top_k)

  top_links <- get_top_links(data, top_k = target_top_k)
  selected_genes <- top_links$link

  print('top selected genes')
  print(head(selected_genes))
  print('number of selected genes:')
  print(length(selected_genes))
 
      
  data[1:5,1:5]
      
  data_filtered  <- data[which(rownames(data) %in% selected_genes),] 
  print('filtered data')
  print(dim(data_filtered))

  if (nrow(data_filtered) < target_top_k) {
    set.seed(123)
    
    n_extra <- target_top_k - nrow(data_filtered)
  
    filtered_indices <- which(rownames(data_all) %in% rownames(data_filtered))
    
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
      data_all[sampled_extra, , drop = FALSE]
    )
  }
  
  
  
  PPP_binary <- data_filtered
  print('final number')
  print(dim(PPP_binary))
  
  getwd()
  
  filetosave <- paste0('PPPbi', fileNum, '.csv')
  
  setwd(netDir)
  
  print(PPP_binary[1:4,1:4])
  PPP_binary <- round(PPP_binary, 4)
  print(PPP_binary[1:4,1:4])
  write.csv(PPP_binary, file = filetosave)
  
  result_dims <- rbind(
    result_dims,
    data.frame(
      matrix_file = file,
      nrow = nrow(PPP_binary),
      ncol = ncol(PPP_binary),
      stringsAsFactors = FALSE
    )
  )

  
}

print("summary of saved PPP_binary dimensions")
print(result_dims)
