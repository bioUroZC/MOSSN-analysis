rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


suppressPackageStartupMessages({
  library(igraph)
})

base_dir <- paste0(PROJ_ROOT, "/1NT/7casestudy")
out_dir <- file.path(base_dir, "output")

edge_file <- file.path(out_dir, "case_display_edges.csv")
node_file <- file.path(out_dir, "case_display_nodes.csv")
summary_file <- file.path(out_dir, "case_selection_summary.csv")

if (!file.exists(edge_file) || !file.exists(node_file) || !file.exists(summary_file)) {
  stop("Missing case-study input files. Please run prepare_case.R first.")
}

edges <- read.csv(edge_file, stringsAsFactors = FALSE)
nodes <- read.csv(node_file, stringsAsFactors = FALSE)
case_summary <- read.csv(summary_file, stringsAsFactors = FALSE)

if (nrow(edges) == 0 || nrow(nodes) == 0) {
  stop("No display edges or nodes were found.")
}

vertex_df <- data.frame(
  name = nodes$gene,
  rewiring_degree = nodes$rewiring_degree,
  rewiring_burden = nodes$rewiring_burden,
  gain_edges = nodes$gain_edges,
  loss_edges = nodes$loss_edges,
  label = nodes$label,
  stringsAsFactors = FALSE
)

all_graph <- graph_from_data_frame(
  d = edges[, c("protein1", "protein2")],
  directed = FALSE,
  vertices = vertex_df
)

layout_xy <- layout_with_fr(all_graph)
rownames(layout_xy) <- V(all_graph)$name

node_scale <- sqrt(vertex_df$rewiring_burden)
node_scale <- 8 + 18 * (node_scale - min(node_scale)) / max(1e-8, max(node_scale) - min(node_scale))
names(node_scale) <- vertex_df$name

node_color_map <- ifelse(vertex_df$gain_edges >= vertex_df$loss_edges, "#c0392b", "#2c7fb8")
names(node_color_map) <- vertex_df$name

edge_width_from_values <- function(values) {
  1.5 + 5 * (values - min(values)) / max(1e-8, max(values) - min(values))
}

plot_panel <- function(panel_edges, edge_weight_col, edge_color, main_title, layout_xy, vertex_df) {
  edge_df <- panel_edges[, c("protein1", "protein2", setdiff(colnames(panel_edges), c("protein1", "protein2"))), drop = FALSE]
  graph_obj <- graph_from_data_frame(
    d = edge_df,
    directed = FALSE,
    vertices = vertex_df
  )
  graph_obj <- delete_vertices(graph_obj, V(graph_obj)[degree(graph_obj) == 0])
  coords <- layout_xy[V(graph_obj)$name, , drop = FALSE]

  edge_values <- edge_attr(graph_obj, edge_weight_col)
  edge_values[is.na(edge_values)] <- 0

  plot(
    graph_obj,
    layout = coords,
    vertex.size = node_scale[V(graph_obj)$name],
    vertex.color = node_color_map[V(graph_obj)$name],
    vertex.frame.color = "white",
    vertex.label = ifelse(V(graph_obj)$label == "", NA, V(graph_obj)$label),
    vertex.label.cex = 0.8,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    edge.width = edge_width_from_values(edge_values),
    edge.color = edge_color,
    main = main_title,
    margin = c(0.2, 0.2, 1, 0.2)
  )
}

normal_edges <- edges[edges$normal_weight > 0, , drop = FALSE]
tumor_edges <- edges[edges$tumor_weight > 0, , drop = FALSE]
differential_edges <- edges

normal_edges$panel_weight <- normal_edges$normal_weight
tumor_edges$panel_weight <- tumor_edges$tumor_weight
differential_edges$panel_weight <- differential_edges$abs_delta
differential_edges$panel_color <- ifelse(differential_edges$direction == "gain", "#d7301f", "#3182bd")

pdf(file.path(out_dir, "case_network_triptych.pdf"), width = 16, height = 6)
par(mfrow = c(1, 3), mar = c(0.5, 0.5, 2.5, 0.5))

plot_panel(
  panel_edges = normal_edges,
  edge_weight_col = "panel_weight",
  edge_color = "#969696",
  main_title = "Matched Normal Network",
  layout_xy = layout_xy,
  vertex_df = vertex_df
)

plot_panel(
  panel_edges = tumor_edges,
  edge_weight_col = "panel_weight",
  edge_color = "#525252",
  main_title = "Matched Tumor Network",
  layout_xy = layout_xy,
  vertex_df = vertex_df
)

diff_graph <- graph_from_data_frame(
  d = differential_edges[, c("protein1", "protein2", setdiff(colnames(differential_edges), c("protein1", "protein2"))), drop = FALSE],
  directed = FALSE,
  vertices = vertex_df
)
coords <- layout_xy[V(diff_graph)$name, , drop = FALSE]

plot(
  diff_graph,
  layout = coords,
  vertex.size = node_scale[V(diff_graph)$name],
  vertex.color = node_color_map[V(diff_graph)$name],
  vertex.frame.color = "white",
  vertex.label = ifelse(V(diff_graph)$label == "", NA, V(diff_graph)$label),
  vertex.label.cex = 0.8,
  vertex.label.family = "sans",
  vertex.label.color = "black",
  edge.width = edge_width_from_values(edge_attr(diff_graph, "abs_delta")),
  edge.color = ifelse(edge_attr(diff_graph, "direction") == "gain", "#d7301f", "#3182bd"),
  main = "Differential Rewiring Network",
  margin = c(0.2, 0.2, 1, 0.2)
)

legend(
  "topleft",
  legend = c("Tumor-gained", "Tumor-lost"),
  col = c("#d7301f", "#3182bd"),
  lwd = 3,
  bty = "n",
  cex = 0.85
)

dev.off()

case_summary_table <- data.frame(
  item = c(
    "Patient",
    "Tumor sample",
    "Normal sample",
    "Distance rank",
    "Displayed nodes",
    "Displayed edges",
    "Tumor-gained edges",
    "Tumor-lost edges",
    "Top rewiring hubs"
  ),
  value = c(
    case_summary$patient[1],
    case_summary$tumor_sample[1],
    case_summary$normal_sample[1],
    paste0(case_summary$edge_distance_rank[1], "/", case_summary$total_patients[1]),
    case_summary$displayed_nodes[1],
    case_summary$displayed_edges[1],
    case_summary$gain_edges_in_display[1],
    case_summary$loss_edges_in_display[1],
    paste(head(nodes$gene[order(nodes$rewiring_burden, decreasing = TRUE)], 5), collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

write.csv(case_summary_table, file.path(out_dir, "case_summary_table.csv"), row.names = FALSE)

cat("Saved case-study figure to:", file.path(out_dir, "case_network_triptych.pdf"), "\n")
cat("Saved summary table to:", file.path(out_dir, "case_summary_table.csv"), "\n")
