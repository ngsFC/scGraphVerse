library(tidyverse)
setwd("/home/francescoc/Desktop/scGRN_simulation/data")

data <- read_tsv("Tcell_string_interactions.tsv")
genes <- unique(c(data$`#node1`, data$node2))
pxp_matrix <- matrix(0, nrow = length(genes), ncol = length(genes))
rownames(pxp_matrix) <- genes
colnames(pxp_matrix) <- genes

for (i in 1:nrow(data)) {
  gene1 <- data$`#node1`[i]
  gene2 <- data$node2[i]
  score <- data$combined_score[i]
  
  
  pxp_matrix[gene1, gene2] <- score
  pxp_matrix[gene2, gene1] <- score 
}

diag(pxp_matrix) <- 1
write.csv(pxp_matrix, "Tcell_weighted_adjacency_matrix.csv", row.names = TRUE)

binary_matrix <- pxp_matrix
binary_matrix[binary_matrix != 0] <- 1

write.csv(binary_matrix, "Tcell_adjacency_matrix.csv", row.names = TRUE)


data <- read_tsv("allblood_string_interactions.tsv")
genes <- unique(c(data$`#node1`, data$node2))
pxp_matrix <- matrix(0, nrow = length(genes), ncol = length(genes))
rownames(pxp_matrix) <- genes
colnames(pxp_matrix) <- genes

for (i in 1:nrow(data)) {
  gene1 <- data$`#node1`[i]
  gene2 <- data$node2[i]
  score <- data$combined_score[i]
  
  pxp_matrix[gene1, gene2] <- score
  pxp_matrix[gene2, gene1] <- score
}

diag(pxp_matrix) <- 1
write.csv(pxp_matrix, "allblood_weighted_adjacency_matrix.csv", row.names = TRUE)

binary_matrix <- pxp_matrix
binary_matrix[binary_matrix != 0] <- 1

write.csv(binary_matrix, "allblood_adjacency_matrix.csv", row.names = TRUE)

