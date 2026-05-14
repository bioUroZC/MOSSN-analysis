# =================================================================

# =================================================================

rm(list = ls())
set.seed(42)


dieasename <- "Paclitaxel"
dataset_name <- "GSE66305"
organ <- "Breast"

base_path <- file.path("/proj/c.zihao/work1/2drugs/",
                       dieasename, "/", 
                       dataset_name, "data")


save_path <- file.path("/proj/c.zihao/work1/2drugs/",
                       dieasename, "/", 
                       dataset_name, "SSP")

dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(save_path, full.names = TRUE, recursive = FALSE), recursive = TRUE, force = TRUE)


setwd(base_path)
cancer <- read.csv("exprSet_filtered.csv", header = TRUE, row.names = 1)
cancer <- cancer[apply(cancer, 1, sd) > 0, ]
colnames(cancer) <- paste0("Cancer_", colnames(cancer))


Normal <- read.csv("/proj/c.zihao/work1/0ref/GTEx/combined_expr_df.csv", row.names = 1)
table(Normal$organ)
Normal <- Normal[which(Normal$organ == organ), ]
Normal$organ <- NULL
Normal <- as.data.frame(t(Normal))
Normal <- Normal[apply(Normal, 1, sd) > 0, ]
colnames(Normal) <- paste0("Normal_", 1:ncol(Normal))


cancer$genes <- rownames(cancer)
Normal$genes <- rownames(Normal)
expr_data <- merge(Normal, cancer, by = "genes")
rownames(expr_data) <- expr_data$genes
expr_data$genes <- NULL


ppi <- read.csv("/proj/c.zihao/work1/1survival/links.csv", header = TRUE, row.names = 1)
head(ppi)
ppi <- ppi[ppi$protein1 %in% rownames(expr_data) & ppi$protein2 %in% rownames(expr_data), ]
network <- ppi[, c("protein1", "protein2")]
colnames(network) <- c("gene1", "gene2")
head(ppi)

# =================================================================

# =================================================================


rank.matrix <- function(x){
  rankmatrix <- sapply(1:ncol(x), function(i) rank(x[, i]))
  colnames(rankmatrix) <- colnames(x)
  rownames(rankmatrix) <- rownames(x)
  return(rankmatrix)
}
ranked_data <- rank.matrix(expr_data)

n_normal <- sum(grepl("Normal_", colnames(ranked_data)))
n_cancer <- sum(grepl("Cancer_", colnames(ranked_data)))

ranked_normal <- ranked_data[, 1:n_normal]
ranked_cancer <- ranked_data[, (n_normal+1):(n_normal+n_cancer)]
mean_rank <- apply(ranked_normal, 1, mean)
combined_rank <- cbind(ranked_normal, mean_rank, ranked_cancer)

# =================================================================

# =================================================================


is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

delta.rank <- function(net, x, n_normal, n_cancer){
  deltarank <- lapply(1:nrow(net), function(i){
    r1 = which(rownames(x) == net[i,1])
    r2 = which(rownames(x) == net[i,2])
    if((length(r1) != 0) & (length(r2) != 0)){
      r = x[r1,] - x[r2,]
      return(c(as.character(net[i,1]), as.character(net[i,2]), r))
    } else {
      return(NULL)
    }
  })
  
  deltarank <- rmNullObs(deltarank)
  deltarank <- do.call(rbind, deltarank)
  
  net.edge = matrix(unlist(deltarank[, 1:2]), ncol=2)
  net.data = as.numeric(deltarank[, n_normal + 3])
  deltarank_cancer = matrix(
    as.numeric(unlist(deltarank[, (n_normal + 4):ncol(deltarank)])),
    ncol = n_cancer
  )
  colnames(deltarank_cancer) <- colnames(x)[(n_normal + 2):(n_normal + 1 + n_cancer)]
  
  return(list(net.edge = net.edge,
              net.data = net.data,
              deltarank_cancer = deltarank_cancer))
}

dresult <- delta.rank(network, combined_rank, n_normal, n_cancer)

# =================================================================

# =================================================================

EPm <- function(net.data, deltarank){
  EPm <- sapply(1:ncol(deltarank), function(i){
    delta <- deltarank[, i] - net.data
    return(delta)
  })
  colnames(EPm) <- colnames(deltarank)
  return(EPm)
}
epm <- EPm(dresult$net.data, dresult$deltarank_cancer)


edge_names <- apply(dresult$net.edge, 1, function(x) paste0(x[1], "_", x[2]))
epm_df <- as.data.frame(epm)
rownames(epm_df) <- edge_names

print(epm_df[1:9,1:9])

colnames(epm_df) <- gsub("^Cancer_", "", colnames(epm_df))

colnames(epm_df)

setwd(save_path)

write.csv(epm_df, "EPm_matrix.csv")

print("All sample for SSP completed")

