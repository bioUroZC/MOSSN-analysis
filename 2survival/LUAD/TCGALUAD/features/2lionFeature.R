# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(igraph)
library(e1071)



projetcNa <- "TCGALUAD"
filenameR <- paste0(paste0(PROJ_ROOT, "/1survival/LUAD/"), projetcNa, "/LION")
filenameR
setwd(filenameR)


df <- read.csv("result.csv", header = T)
df[1:5,1:5]


df$Node1_clean <- pmin(df$reg, df$tar)
df$Node2_clean <- pmax(df$reg, df$tar)
df$link <- paste0(df$Node1_clean, "_", df$Node2_clean)

df[1:5,1:5]

df$reg <- NULL
df$tar <- NULL
df$Node1_clean <- NULL
df$Node2_clean <- NULL

rownames(df) <- df$link

df$link  <- NULL

df[df<0.7] <- 0

df$link <- rownames(df)

df <- df %>%
  dplyr::select(link, everything())

df[1:5,1:5]


# ===================================================

# ===================================================

PPI <- df


# Initialize results data frame
results <- data.frame()

# Extract Node1 and Node2 from 'link' column
edges <- do.call(rbind, strsplit(PPI$link, "_"))
colnames(edges) <- c("Node1", "Node2")

# Loop over each sample (column 2 onward)
for (i in 2:ncol(PPI)) {
  sample_name <- colnames(PPI)[i]
  print(paste("Processing:", sample_name, "(Column", i, ")"))
  sample_weights <- PPI[[i]]
  
  # Create edge list with current sample weights
  data <- data.frame(
    Node1 = edges[, 1],
    Node2 = edges[, 2],
    UpdatedWeight = sample_weights
  )
  
  # Filter for edges with weight > 0.7
  data <- subset(data, UpdatedWeight > 0.7)
  data$UpdatedWeight <- round(data$UpdatedWeight, 3)
  
  # Skip if too few edges
  if (nrow(data) < 2) next
  
  # Create graph
  g <- graph_from_data_frame(data, directed = FALSE)
  E(g)$weight <- data$UpdatedWeight
  degrees <- degree(g)
  
  # Network-level features
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
  
  # Edge weight features
  avg_edge_weight <- mean(data$UpdatedWeight)
  max_edge_weight <- max(data$UpdatedWeight)
  min_edge_weight <- min(data$UpdatedWeight)
  edge_weight_sd <- sd(data$UpdatedWeight)
  
  # Centralities
  pagerank <- page_rank(g, weights = E(g)$weight)$vector
  mean_pagerank <- mean(pagerank)
  sd_pagerank <- sd(pagerank)
  max_pagerank <- max(pagerank)
  
  eigen_centrality <- eigen_centrality(g, directed = FALSE, weights = E(g)$weight)$vector
  mean_eigen_centrality <- mean(eigen_centrality)
  sd_eigen_centrality <- sd(eigen_centrality)
  max_eigen_centrality <- max(eigen_centrality)
  
  # Largest component
  comp <- components(g)
  largest_component_nodes <- which(comp$membership == which.max(comp$csize))
  largest_component_subgraph <- induced_subgraph(g, largest_component_nodes)
  
  largest_component_size <- length(largest_component_nodes)
  largest_component_edges <- ecount(largest_component_subgraph)
  largest_component_density <- edge_density(largest_component_subgraph)
  largest_component_clustering <- transitivity(largest_component_subgraph, type = "global")
  connected_components <- comp$no
  
  # Sample top PageRank nodes (up to 1000) for path analysis
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
  
  # Combine all features into one row
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

# Final result
print(results[1:6,1:6])



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

filenameR <- paste0(paste0(PROJ_ROOT, "/1survival/LUAD/"), projetcNa, "/features")
filenameR
setwd(filenameR)

write.csv(wrwData, file = 'lionData.csv')
cat("files saved")

