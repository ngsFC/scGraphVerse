# NodeVerse: A Package for Network Inference

## Overview
NodeVerse is an R package designed for inferring and analyzing gene regulatory networks (GRNs). The package provides tools for generating adjacency matrices, inferring networks using different methods, comparing inferred networks to ground truth, and visualizing network properties.

## Installation
To install NodeVerse, you can use the following command:

```r
# Install from GitHub (if hosted there)
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("ngsFC/NodeVerse")
```

## Functions

### 1. `infer_networks()`
Infers gene regulatory networks using different methods, including GENIE3, GRNBoost2, ZILGM, JRF, and PCzinb.

#### Usage:
```r
network <- infer_networks(count_matrices_list, method = "GENIE3")
```

### 2. `generate_adjacency()`
Generates adjacency matrices from a list of data frames containing gene interaction data.

#### Usage:
```r
adjacency_matrices <- generate_adjacency(data_frames_list)
```

### 3. `symmetrize()`
Ensures adjacency matrices are symmetric by applying a specified function to corresponding off-diagonal elements.

#### Usage:
```r
symmetric_matrices <- symmetrize(matrix_list, method = "mean")
```

### 4. `plotROC()`
Plots the Receiver Operating Characteristic (ROC) curve for different inferred networks.

#### Usage:
```r
plotROC(predicted_matrices, ground_truth)
```

### 5. `cutoff_adjacency()`
Applies a cutoff to adjacency matrices based on percentile values from shuffled networks.

#### Usage:
```r
filtered_adjacency <- cutoff_adjacency(count_matrices, method = "GRNBoost2")
```

### 6. `pscores()`
Computes performance metrics for predicted adjacency matrices compared to a ground truth.

#### Usage:
```r
pscores(predicted_matrices, ground_truth)
```

### 7. `plotg()`
Generates network visualizations from adjacency matrices.

#### Usage:
```r
plotg(adjacency_matrices)
```

### 8. `create_consensus()`
Creates a consensus adjacency matrix using voting, union, or INet methods.

#### Usage:
```r
consensus_matrix <- create_consensus(adjacency_list, method = "vote")
```

### 9. `compare_consensus()`
Compares consensus adjacency matrices to original graphs and visualizes their differences.

#### Usage:
```r
compare_consensus(ground_truth_matrix, consensus_matrix)
```

## License
This package is licensed under the MIT License.

## Aknowledge
"National Centre for HPC, Big Data and Quantum Computing" - CN00000013 - CUP B93C22999629996

## Author
Developed by [Francesco Cecere].

