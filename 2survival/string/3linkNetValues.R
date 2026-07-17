# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)
library(igraph)
library(poweRlaw)

# ===================================================

# ===================================================

setwd(paste0(PROJ_ROOT, "/1survival"))

ppi_links <- read.csv("links.csv", row.names = 1)
ppi_links <- subset(ppi_links, score >= 0.9)  


set.seed(123)  # 保证可复现
n_edges <- nrow(ppi_links)
n_sample <- ceiling(n_edges * 0.10)
random_links <- ppi_links[sample(n_edges, n_sample), ]

ppi_links <- random_links 

# 2. 构建无向图
g <- graph_from_data_frame(ppi_links[, c("protein1", "protein2")], directed = FALSE)

# 3. 幂律指数 gamma
deg <- degree(g)
gamma <- NA
if(length(unique(deg)) > 1){
  pl_fit <- try(displ$new(deg + 1), silent = TRUE)
  if(!inherits(pl_fit, "try-error")){
    est <- estimate_xmin(pl_fit)
    pl_fit$setXmin(est)
    gamma <- pl_fit$pars
  }
}
cat("Power-law gamma:", gamma, "\n")

# 4. 聚类系数 clustering coefficient
clustering <- transitivity(g, type = "average")
cat("Clustering coefficient:", clustering, "\n")

# 5. 平均最短路径 average shortest path
avg_path <- NA
comp <- components(g)
if(comp$no > 0 && max(comp$csize) >= 2){
  largest_comp <- which(comp$membership == which.max(comp$csize))
  subg <- induced_subgraph(g, largest_comp)
  avg_path <- average.path.length(subg)
}
cat("Average shortest path:", avg_path, "\n")

# 6. Louvain模块度 modularity
mod <- NA
if(ecount(g) > 1){
  comm <- try(cluster_louvain(g), silent = TRUE)
  if(!inherits(comm, "try-error")){
    mod <- modularity(comm)
  }
}
cat("Louvain modularity:", mod, "\n")
