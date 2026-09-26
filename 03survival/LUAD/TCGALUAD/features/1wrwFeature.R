rm(list = ls())

library(dplyr)
library(igraph)
library(e1071)

projetcNa <- "TCGALUAD"
filenameR <- paste0("/proj/c.zihao/work1/03survival/LUAD/", projetcNa, "/WRW")
filenameR
setwd(filenameR)

edges_files <- list.files(path = ".", pattern = "edges\\.csv$", recursive = TRUE, full.names = TRUE)

merged_data <- NULL
total_files <- length(edges_files)

for (i in seq_along(edges_files)) {
  file <- edges_files[i]
  
  cat(sprintf("[%d/%d] (%.1f%%) Processing: %s\n", i, total_files, i / total_files * 100, file))
  
  df <- read.csv(file)
  df <- subset(df, select = c("Node1", "Node2", "UpdatedWeight"))
  df <- subset(df, UpdatedWeight > 0.7)
  df$Node1_clean <- pmin(df$Node1, df$Node2)
  df$Node2_clean <- pmax(df$Node1, df$Node2)
  df$link <- paste0(df$Node1_clean, "_", df$Node2_clean)
  
  file_name <- gsub("^\\./|_edges\\.csv$", "", file)
  df <- df[, c("link", "UpdatedWeight")]
  colnames(df)[2] <- file_name
  
  if (is.null(merged_data)) {
    merged_data <- df
  } else {
    merged_data <- merge(merged_data, df, by = "link", all = TRUE)
  }
}

head(merged_data)[1:5,1:5]

merged_data[is.na(merged_data)] <- 0

PPI <- merged_data

results <- data.frame()

edges <- do.call(rbind, strsplit(PPI$link, "_"))
colnames(edges) <- c("Node1", "Node2")

for (i in 2:ncol(PPI)) {
  sample_name <- colnames(PPI)[i]
  print(paste("Processing:", sample_name, "(Column", i, ")"))
  sample_weights <- PPI[[i]]
  
  data <- data.frame(
    Node1 = edges[, 1],
    Node2 = edges[, 2],
    UpdatedWeight = sample_weights
  )
  
  data <- subset(data, UpdatedWeight > 0.7)
  data$UpdatedWeight <- round(data$UpdatedWeight, 3)
  
  if (nrow(data) < 2) next
  
  g <- graph_from_data_frame(data, directed = FALSE)
  E(g)$weight <- data$UpdatedWeight
  degrees <- degree(g)
  
  num_nodes <- vcount(g)
  num_edges <- ecount(g)
  avg_degree <- mean(degrees)
  std_degree <- sd(degrees)
  max_degree <- max(degrees)
  min_degree <- min(degrees)
  graph_density <- edge_density(g)
  clustering_coeff_global <- transitivity(g, type = "global")
  clustering_coeff_avg <- mean(transitivity(g, type = "local"), na.rm = TRUE)
  avg_neighbor_degree <- mean(knn(g)$knn, na.rm = TRUE)
  degree_heterogeneity <- var(degrees) / avg_degree^2
  degree_skewness <- skewness(degrees)
  degree_kurtosis <- kurtosis(degrees)
  assortativity <- assortativity_degree(g)
  
  avg_edge_weight <- mean(data$UpdatedWeight)
  max_edge_weight <- max(data$UpdatedWeight)
  min_edge_weight <- min(data$UpdatedWeight)
  edge_weight_sd <- sd(data$UpdatedWeight)
  
  pagerank <- page_rank(g, weights = E(g)$weight)$vector
  mean_pagerank <- mean(pagerank)
  sd_pagerank <- sd(pagerank)
  max_pagerank <- max(pagerank)
  
  eigen_centrality <- eigen_centrality(g, directed = FALSE, weights = E(g)$weight)$vector
  mean_eigen_centrality <- mean(eigen_centrality)
  sd_eigen_centrality <- sd(eigen_centrality)
  max_eigen_centrality <- max(eigen_centrality)
  
  comp <- components(g)
  largest_component_nodes <- which(comp$membership == which.max(comp$csize))
  largest_component_subgraph <- induced_subgraph(g, largest_component_nodes)
  
  largest_component_size <- length(largest_component_nodes)
  largest_component_edges <- ecount(largest_component_subgraph)
  largest_component_density <- edge_density(largest_component_subgraph)
  largest_component_clustering <- transitivity(largest_component_subgraph, type = "global")
  connected_components <- comp$no
  
  sub_pr <- page_rank(largest_component_subgraph, weights = E(largest_component_subgraph)$weight)$vector
  top_nodes <- order(sub_pr, decreasing = TRUE)[1:min(1000, vcount(largest_component_subgraph))]
  sampled_subgraph <- induced_subgraph(largest_component_subgraph, top_nodes)
  
  distances_mat <- distances(sampled_subgraph, weights = E(sampled_subgraph)$weight)
  finite_distances <- distances_mat[distances_mat > 0 & is.finite(distances_mat)]
  
  diameter <- if (length(finite_distances) > 0) max(finite_distances) else NA
  avg_path_length <- if (length(finite_distances) > 0) mean(finite_distances) else NA
  global_efficiency <- if (length(finite_distances) > 0) mean(1 / finite_distances) else NA
  
  eccentricity_values <- apply(distances_mat, 1, max, na.rm = TRUE)
  eccentricity_values <- eccentricity_values[is.finite(eccentricity_values)]
  eccentricity_avg <- if (length(eccentricity_values) > 0) mean(eccentricity_values) else NA
  radius <- if (length(eccentricity_values) > 0) min(eccentricity_values) else NA
  
  sample_features <- data.frame(
    Sample = sample_name,
    NumNodes = num_nodes,
    NumEdges = num_edges,
    AvgDegree = avg_degree,
    StdDegree = std_degree,
    MaxDegree = max_degree,
    MinDegree = min_degree,
    Density = graph_density,
    GlobalClustering = clustering_coeff_global,
    AvgLocalClustering = clustering_coeff_avg,
    AvgNeighborDegree = avg_neighbor_degree,
    DegreeHeterogeneity = degree_heterogeneity,
    DegreeSkewness = degree_skewness,
    DegreeKurtosis = degree_kurtosis,
    Assortativity = assortativity,
    AvgEdgeWeight = avg_edge_weight,
    MaxEdgeWeight = max_edge_weight,
    MinEdgeWeight = min_edge_weight,
    EdgeWeightSD = edge_weight_sd,
    MeanPageRank = mean_pagerank,
    SDPageRank = sd_pagerank,
    MaxPageRank = max_pagerank,
    MeanEigenCentrality = mean_eigen_centrality,
    SDEigenCentrality = sd_eigen_centrality,
    MaxEigenCentrality = max_eigen_centrality,
    LargestComponentSize = largest_component_size,
    LargestComponentEdges = largest_component_edges,
    LargestComponentDensity = largest_component_density,
    LargestComponentClustering = largest_component_clustering,
    Diameter = diameter,
    AvgPathLength = avg_path_length,
    GlobalEfficiency = global_efficiency,
    EccentricityAvg = eccentricity_avg,
    Radius = radius,
    ConnectedComponents = connected_components
  )
  
  results <- rbind(results, sample_features)
}

print(results[1:6,1:6])

filenameR <- paste0("/proj/c.zihao/work1/03survival/LUAD/", projetcNa, "/features")
filenameR
setwd(filenameR)

write.csv(results, file = 'Netfeature.csv')
cat("files saved")

PPI$rowmean <- rowMeans(PPI[ , -1])
PPI_top <- PPI[order(-PPI$rowmean), ][1:500, ]
PPI_top$rowmean <- NULL
rownames(PPI_top) <- PPI_top$link
PPI_top$link <- NULL
PPI_top <- as.data.frame(t(PPI_top))
head(PPI_top)[1:6,1:6]
PPI_top$Sample <- rownames(PPI_top)

print(PPI_top[1:6,1:6])

wrwData <- merge(results, PPI_top, by="Sample")
print(wrwData[1:6,1:6])

filenameR <- paste0("/proj/c.zihao/work1/03survival/LUAD/", projetcNa, "/features")
filenameR
setwd(filenameR)
write.csv(wrwData, file = 'wrwData.csv')
cat("files saved")
