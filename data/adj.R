# Load necessary library
library(dplyr)
setwd("/home/francescoc/Desktop/scGRN_simulation/data")

# Read the data
file_path <- "string_interactions.tsv"  # Replace with your actual file path
data <- read.table(file_path, header = F, sep = "\t", stringsAsFactors = FALSE)

# Assign column names to the data frame
colnames(data) <- c("Gene1","Gene2","node1_string_id","node2_string_id","neighborhood_on_chromosome","gene_fusion","phylogenetic_cooccurrence","homology","coexpression","experimentally_determined_interaction","database_annotated","automated_textmining","combined_score")

# Create a list of unique genes
genes <- unique(c(data$Gene1, data$Gene2))

# Initialize an adjacency matrix with zeros
adj_matrix <- matrix(0, nrow = length(genes), ncol = length(genes), 
                     dimnames = list(genes, genes))

# Populate the adjacency matrix with weights
for (i in 1:nrow(data)) {
  gene1 <- data$Gene1[i]
  gene2 <- data$Gene2[i]
  weight <- data$combined_score[i]
  
  # Assign weight to both gene pairs (symmetrically)
  adj_matrix[gene1, gene2] <- weight
  adj_matrix[gene2, gene1] <- weight
}

# Convert the adjacency matrix to a data frame for easier viewing
adj_matrix_df <- as.data.frame(adj_matrix)

# Print the adjacency matrix
print(adj_matrix_df)

# Save the adjacency matrix to a CSV file (optional)
write.csv(adj_matrix_df, "weighted_adjacency_matrix.csv", row.names = TRUE)

# Convert all values greater than 0 to 1
adj_matrix_df[adj_matrix_df > 0] <- 1

# Save the transformed data to a new CSV file
write.csv(adj_matrix_df, "adjacency_matrix.csv", row.names = TRUE)
