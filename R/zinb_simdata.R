#' Simulate Multiple Zero-Inflated Negative Binomial (ZINB) Count Matrices with Cell-Specific Depth
#'
#' This function generates multiple count matrices following a **zero-inflated negative binomial (ZINB)**
#' distribution using an adjacency matrix that defines gene-gene relationships. Each matrix simulates 
#' gene expression across a specified number of cells, incorporating sequencing depth variability.
#'
#' @param n Integer. The number of samples (cells) in each simulated matrix.
#' @param p Integer. The number of variables (genes) in each matrix.
#' @param B Matrix. A **symmetric binary adjacency matrix** (`0/1`), where rows and columns correspond 
#'   to gene names. The adjacency matrix defines gene-gene connectivity, influencing expression levels.
#' @param mu_range List of numeric vectors. Each element is a numeric vector of length 2, specifying 
#'   the range of mean expression (`mu`) values per matrix. The function draws gene-specific means from 
#'   a uniform distribution within these ranges.
#' @param mu_noise Numeric vector of length `kmat`. Mean of the background noise component for each matrix.
#' @param theta Numeric vector of length `kmat`. The dispersion parameter of the **negative binomial distribution** 
#'   for each matrix. A lower value of `theta` increases overdispersion.
#' @param pi Numeric vector of length `kmat`. The probability of excess zeros (`0 < pi < 1`) for each matrix, 
#'   controlling sparsity due to dropouts or biological factors.
#' @param kmat Integer. The number of count matrices to generate (default: `1`).
#' @param depth_range Numeric vector of length 2 or `NA`. If provided, defines the range of total sequencing 
#'   depth per cell (e.g., `c(500, 5000)`). Default: `NA` (no depth scaling applied).
#' @param seed Integer (optional). A random seed for reproducibility.
#' 
#' @return A **list** of `kmat` numeric matrices, each with dimensions `n x p`, where:
#'   - **Rows** correspond to cells (`cell_1, cell_2, ..., cell_n`).
#'   - **Columns** correspond to genes (from `rownames(B)`).
#'   - Expression values follow a **zero-inflated negative binomial (ZINB)** distribution.
#'
#' @details 
#' The function simulates count data by:
#' 1. Generating gene expression values from a **zero-inflated negative binomial (ZINB)** model.
#' 2. Modifying gene-gene expression levels based on the adjacency matrix `B`.
#' 3. Incorporating **cell-specific sequencing depth** if `depth_range` is specified.
#'
#' **Biological Interpretation:**
#' - This method is useful for benchmarking **single-cell RNA sequencing (scRNA-seq)** algorithms, 
#'   where excess zeros may arise due to **dropout events**.
#' - The adjacency matrix `B` allows simulation of **gene regulatory network (GRN)** interactions.
#' - Sequencing depth variation mimics experimental conditions in **scRNA-seq datasets**.
#'
#' @importFrom distributions3 rzinbinom
#' @examples
#' \dontrun{
#' # Define adjacency matrix (example with 5 genes)
#' B <- matrix(0, nrow = 5, ncol = 5)
#' rownames(B) <- colnames(B) <- paste0("Gene", 1:5)
#' B[1, 2] <- B[2, 3] <- B[4, 5] <- 1  # Define some gene interactions
#' B <- B + t(B)  # Make symmetric
#'
#' # Simulate data with 100 cells, 5 genes, and 2 matrices
#' zinb_matrices <- zinb_simdata(
#'   n = 100, p = 5, B = B,
#'   mu_range = list(c(1, 5), c(2, 6)), mu_noise = c(0.5, 0.7),
#'   theta = c(1, 2), pi = c(0.2, 0.3), kmat = 2, 
#'   depth_range = c(500, 5000), seed = 42
#' )
#'
#' # View one of the simulated matrices
#' print(zinb_matrices[[1]])
#' }
#' @export
zinb_simdata <- function(n, p, B, mu_range, mu_noise, theta, pi, kmat = 1, depth_range = NA, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(is.numeric(n), n > 0, floor(n) == n)
  stopifnot(is.numeric(p), p > 0, floor(p) == p)
  stopifnot(is.matrix(B), nrow(B) == ncol(B))
  stopifnot(all(B %in% c(0, 1)))
  stopifnot(is.numeric(kmat), kmat > 0, floor(kmat) == kmat)
  stopifnot(length(mu_range) == kmat, all(sapply(mu_range, function(x) length(x) == 2 && all(x > 0))))
  stopifnot(length(mu_noise) == kmat, all(mu_noise >= 0))
  stopifnot(length(theta) == kmat, all(theta > 0))
  stopifnot(length(pi) == kmat, all(pi > 0 & pi < 1))
  
  if (!is.na(depth_range[1])) {
    stopifnot(is.numeric(depth_range), length(depth_range) == 2, all(depth_range > 0), depth_range[1] < depth_range[2])
  }

  gene_names <- rownames(B)
  cellID <- paste0("cell_", seq_len(n))
  
  B <- ifelse(B > 0, 1, 0)
  
  edges <- which(B == 1, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], ] 
  
  A <- diag(1, nrow = p, ncol = p)
  for (i in seq_len(nrow(edges))) {
    tmp <- rep(0, p)
    tmp[edges[i, ]] <- 1
    A <- cbind(A, tmp)
  }
  
  B[edges] <- sample(seq(min(unlist(mu_range)), max(unlist(mu_range))), length(edges[, 1]), replace = TRUE)
  B <- (B | t(B)) * 1 
  
  matrices <- vector("list", kmat)
  
  for (k in seq_len(kmat)) {
    mu <- runif(p, mu_range[[k]][1], mu_range[[k]][2])
    
    sigma <- B  
    nonzero_sigma <- sigma[lower.tri(sigma) & sigma != 0]
    Y_mu <- c(mu, nonzero_sigma)  
    
    Y <- matrix(rzinbinom(length(Y_mu) * n, mu = rep(Y_mu, each = n), theta = theta[k], pi = pi[k]), 
                nrow = length(Y_mu), ncol = n)
    X <- A %*% Y
    
    noise_matrix <- matrix(rzinbinom(n * p, mu = mu_noise[k], theta = 1, pi = pi[k]), nrow = p, ncol = n)
    X <- X + noise_matrix
    
    X <- t(X)
    if (!is.null(gene_names)) colnames(X) <- gene_names
    rownames(X) <- cellID
    
    if (!is.na(depth_range[1])) {
      cell_depths <- runif(n, min = depth_range[1], max = depth_range[2])
      
      row_sums <- rowSums(X)
      row_sums[row_sums == 0] <- 1
      
      X <- sweep(X, 1, row_sums, FUN = "/")
      X[is.na(X)] <- 0
      
      X <- sweep(X, 1, cell_depths, FUN = "*")
      X <- round(X)
    }
    
    matrices[[k]] <- X
  }
  
  return(matrices)
}

