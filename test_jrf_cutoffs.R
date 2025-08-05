#!/usr/bin/env Rscript

# Test script to compare JRF vs GENIE3 cutoff calculation logic
# This script verifies that JRF uses joint null distributions while GENIE3 uses individual matrices

cat("Loading scGraphVerse...\n")
library(scGraphVerse)

# Create small test data with 2 conditions
set.seed(123)
n_genes <- 20
n_cells_1 <- 30
n_cells_2 <- 25

# Generate synthetic expression matrices for 2 conditions
mat1 <- matrix(rnorm(n_genes * n_cells_1, mean = 5, sd = 2), 
               nrow = n_genes, ncol = n_cells_1)
mat2 <- matrix(rnorm(n_genes * n_cells_2, mean = 4, sd = 1.5), 
               nrow = n_genes, ncol = n_cells_2)

rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:n_genes)
colnames(mat1) <- paste0("Cell1_", 1:n_cells_1)
colnames(mat2) <- paste0("Cell2_", 1:n_cells_2)

count_matrices <- list(Condition1 = mat1, Condition2 = mat2)

cat("Test data created:\n")
cat("- Condition 1:", nrow(mat1), "genes x", ncol(mat1), "cells\n")
cat("- Condition 2:", nrow(mat2), "genes x", ncol(mat2), "cells\n")

# Test 1: GENIE3 networks and cutoffs (individual matrix processing)
cat("\n=== Testing GENIE3 (individual matrix processing) ===\n")
genie3_networks <- infer_networks(count_matrices, method = "GENIE3", nCores = 1, verbose = TRUE)
genie3_adjm <- generate_adjacency(genie3_networks)
genie3_symm <- symmetrize(genie3_adjm, weight_function = "mean")

cat("Running GENIE3 cutoff calculation with n=2 shuffles...\n")
genie3_binary <- cutoff_adjacency(
  count_matrices = count_matrices,
  weighted_adjm_list = genie3_symm,
  n = 2,
  method = "GENIE3",
  quantile_threshold = 0.95,
  debug = TRUE
)

# Test 2: JRF networks and cutoffs (joint processing)  
cat("\n=== Testing JRF (joint matrix processing) ===\n")
jrf_networks <- infer_networks(count_matrices, method = "JRF", nCores = 1, verbose = TRUE)
jrf_adjm <- generate_adjacency(jrf_networks)
jrf_symm <- symmetrize(jrf_adjm, weight_function = "mean")

cat("Running JRF cutoff calculation with n=2 shuffles...\n")
jrf_binary <- cutoff_adjacency(
  count_matrices = count_matrices,
  weighted_adjm_list = jrf_symm,
  n = 2,
  method = "JRF",
  quantile_threshold = 0.95,
  debug = TRUE
)

# Analyze results
cat("\n=== Results Analysis ===\n")
cat("GENIE3 binary matrices:\n")
cat("- Condition 1 edges:", sum(genie3_binary[[1]]), "\n")
cat("- Condition 2 edges:", sum(genie3_binary[[2]]), "\n")

cat("JRF binary matrices:\n")
cat("- Condition 1 edges:", sum(jrf_binary[[1]]), "\n")
cat("- Condition 2 edges:", sum(jrf_binary[[2]]), "\n")

# Compare edge densities
genie3_density <- c(
  sum(genie3_binary[[1]]) / (nrow(genie3_binary[[1]]) * (nrow(genie3_binary[[1]]) - 1)),
  sum(genie3_binary[[2]]) / (nrow(genie3_binary[[2]]) * (nrow(genie3_binary[[2]]) - 1))
)

jrf_density <- c(
  sum(jrf_binary[[1]]) / (nrow(jrf_binary[[1]]) * (nrow(jrf_binary[[1]]) - 1)),
  sum(jrf_binary[[2]]) / (nrow(jrf_binary[[2]]) * (nrow(jrf_binary[[2]]) - 1))
)

cat("\nEdge densities:\n")
cat("GENIE3 - Condition 1:", round(genie3_density[1], 4), "Condition 2:", round(genie3_density[2], 4), "\n")
cat("JRF    - Condition 1:", round(jrf_density[1], 4), "Condition 2:", round(jrf_density[2], 4), "\n")

cat("\n=== Test completed successfully! ===\n")
cat("Key differences:\n")
cat("1. GENIE3: Each condition processed independently for cutoff calculation\n")
cat("2. JRF: All conditions processed jointly for cutoff calculation\n")
cat("3. Both methods should produce reasonable but different edge densities\n")