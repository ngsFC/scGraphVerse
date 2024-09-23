# Load necessary library
library(dplyr)
setwd("/home/francescoc/Desktop/scGRN_simulation/data/genemania-interactions.txt")

# Read the data
file_path <- "genemania-interactions.txt"  # Replace with your actual file path
data <- read.table(file_path, header = T, sep = "\t", stringsAsFactors = FALSE)

# Assign column names to the data frame
colnames(data) <- c("Gene1", "Gene2", "Weight", "NetworkGroup", "Network")

# Create a list of unique genes
genes <- unique(c(data$Gene1, data$Gene2))

# Initialize an adjacency matrix with zeros
adj_matrix <- matrix(0, nrow = length(genes), ncol = length(genes), 
                     dimnames = list(genes, genes))

# Populate the adjacency matrix with weights
for (i in 1:nrow(data)) {
  gene1 <- data$Gene1[i]
  gene2 <- data$Gene2[i]
  weight <- data$Weight[i]
  
  # Assign weight to both gene pairs (symmetrically)
  adj_matrix[gene1, gene2] <- weight
  adj_matrix[gene2, gene1] <- weight
}

# Convert the adjacency matrix to a data frame for easier viewing
adj_matrix_df <- as.data.frame(adj_matrix)

# Print the adjacency matrix
print(adj_matrix_df)

# Save the adjacency matrix to a CSV file (optional)
write.csv(adj_matrix_df, "adjacency_matrix.csv", row.names = TRUE)

