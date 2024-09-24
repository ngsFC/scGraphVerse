# Load required libraries
library(GENIE3)
library(doParallel)
library(igraph)
library(tidyverse)

setwd("/home/francescoc/Desktop/scGRN_simulation/Bnet")
count_matrix <- readRDS("./../data/simatx.RDS")
adjm <- read.table("./../data/adjacency_matrix.csv", header = T, row.names = 1, sep = ",")
marker <- read.table("./../data/Tcell.marker.csv", header = T, sep = ",")

count_matrix1 <- as.data.frame(count_matrix[1])
colnames(count_matrix1) <- colnames(adjm)
count_matrix2 <- as.data.frame(count_matrix[2])
colnames(count_matrix2) <- colnames(adjm)
count_matrix3 <- as.data.frame(count_matrix[3])
colnames(count_matrix3) <- colnames(adjm)
count_matrix4 <- as.data.frame(count_matrix[4])
colnames(count_matrix4) <- colnames(adjm)
count_matrix5 <- as.data.frame(count_matrix[5])
colnames(count_matrix5) <- colnames(adjm)

# Run GENIE3 using the count matrix
set.seed(123)
regulatory_network_genie3 <- GENIE3(t(count_matrix1))

# Extract link list (gene regulatory interactions) from GENIE3 results
link_list_genie3 <- getLinkList(regulatory_network_genie3)

# Save the results
write.csv(link_list_genie3, "genie3_network.csv", row.names = FALSE)
#link_list_genie3 <- link_list_genie3 %>% filter(weight >= 0.1) 

# List of unique genes
gene_names <- unique(c(link_list_genie3$regulator, link_list_genie3$target))

# Create an empty adjacency matrix
adj_matrix_genie3 <- matrix(0, nrow = length(gene_names), ncol = length(gene_names))
rownames(adj_matrix_genie3) <- colnames(adj_matrix_genie3) <- gene_names

# Fill the adjacency matrix based on the links from GENIE3 with a weight condition
for (i in 1:nrow(link_list_genie3)) {
  regulator <- link_list_genie3$regulator[i]
  target <- link_list_genie3$target[i]
  weight <- link_list_genie3$weight[i]
  
  # Only set 1 if the weight is >= 0.1
  if (weight >= 0.1) {
    adj_matrix_genie3[regulator, target] <- 1
  }
}

# Export the adjacency matrix if needed
write.csv(adj_matrix_genie3, "genie3_adjacency_matrix.csv")

# Create igraph objects for both networks
graph_genie3 <- graph_from_adjacency_matrix(adj_matrix_genie3, mode = "undirected")
graph_provided <- graph_from_adjacency_matrix(as.matrix(adjm), mode = "undirected")

# Plot both networks side by side
par(mfrow = c(1, 2))  # Side by side plotting

# Plot GENIE3 Network
plot(graph_genie3, main = "GENIE3 Inferred Network", vertex.label.color = "black",
     vertex.size = 10, edge.arrow.size = 0.5, vertex.label.cex = 0.7)

# Plot Provided Network
plot(graph_provided, main = "Provided Network", vertex.label.color = "black",
     vertex.size = 10, edge.arrow.size = 0.5, vertex.label.cex = 0.7)

par(mfrow = c(1, 1))  # Reset plotting layout

# Get the edges for both networks
edges_genie3 <- igraph::as_data_frame(graph_genie3, what = "edges")
edges_provided <- igraph::as_data_frame(graph_provided, what = "edges")

# Convert edges to sets of pairs (for easy comparison)
edges_genie3_set <- paste(edges_genie3$from, edges_genie3$to, sep = "-")
edges_provided_set <- paste(edges_provided$from, edges_provided$to, sep = "-")

# Calculate the Jaccard similarity between the two sets of edges
common_edges <- intersect(edges_genie3_set, edges_provided_set)
union_edges <- union(edges_genie3_set, edges_provided_set)
jaccard_similarity <- length(common_edges) / length(union_edges)

cat("Jaccard Similarity: ", jaccard_similarity, "\n")

