#' Simulate Zero-Inflated Negative Binomial (ZINB) Count Matrices with Sequencing Depth
#'
#' Simulates one or more count matrices following a zero-inflated negative binomial (ZINB) distribution,
#' incorporating gene-gene interaction structures and cell-specific sequencing depth variation.
#'
#' @param n Integer. Number of cells (samples) in each simulated matrix.
#' @param p Integer. Number of genes (features) in each simulated matrix.
#' @param B A symmetric binary adjacency matrix (0/1) defining gene-gene connectivity.
#'   Row and column names correspond to gene names.
#' @param mu_range List of numeric vectors (length 2 each). Range of gene expression means
#'   for each simulated matrix.
#' @param mu_noise Numeric vector. Mean of background noise for each matrix.
#' @param theta Numeric vector. Dispersion parameters of the negative binomial distribution
#'   for each matrix. Smaller \code{theta} implies higher overdispersion.
#' @param pi Numeric vector. Probability of excess zeros (\code{0 < pi < 1}) for each matrix.
#' @param kmat Integer. Number of count matrices to simulate. Default is \code{1}.
#' @param depth_range Numeric vector of length 2 or \code{NA}. Range of total sequencing depth per cell.
#'   If \code{NA}, no depth adjustment is performed.
#'
#' @return
#' A list containing \code{kmat} matrices. Each matrix has:
#' \itemize{
#'   \item Rows representing cells (\code{cell_1, cell_2, ..., cell_n}).
#'   \item Columns representing genes (based on \code{rownames(B)}).
#'   \item Count values following a ZINB distribution.
#' }
#'
#' @details
#' Each simulated matrix:
#' \enumerate{
#'   \item Generates gene expression values based on a ZINB model.
#'   \item Modulates expression using the input adjacency matrix \code{B}.
#'   \item Applies random sequencing depth scaling if \code{depth_range} is provided.
#' }
#'
#' The function is particularly useful for benchmarking single-cell RNA-seq network inference
#' methods where dropout events and network structure are relevant.
#'
#' @importFrom distributions3 rzinbinom
#' @export
#'
#' @examples
#' # Create a simple adjacency matrix
#' B <- matrix(0, nrow = 5, ncol = 5)
#' rownames(B) <- colnames(B) <- paste0("Gene", 1:5)
#' B[1, 2] <- B[2, 3] <- B[4, 5] <- 1
#' B <- B + t(B) # Ensure symmetry
#'
#' # Simulate two count matrices
#' zinb_matrices <- zinb_simdata(
#'   n = 100, p = 5, B = B,
#'   mu_range = list(c(1, 5), c(2, 6)),
#'   mu_noise = c(0.5, 0.7),
#'   theta = c(1, 2),
#'   pi = c(0.2, 0.3),
#'   kmat = 2,
#'   depth_range = c(500, 5000),
#' )
#'
#' # View one simulated matrix
#' print(zinb_matrices[[1]])
zinb_simdata <- function(n, p, B, mu_range, mu_noise, theta, pi, kmat = 1, depth_range = NA) {
  stopifnot(is.numeric(n), n > 0, floor(n) == n)
  stopifnot(is.numeric(p), p > 0, floor(p) == p)
  stopifnot(is.matrix(B), nrow(B) == ncol(B), all(B %in% c(0, 1)))
  stopifnot(is.numeric(kmat), kmat > 0, floor(kmat) == kmat)
  stopifnot(length(mu_range) == kmat, all(vapply(mu_range, function(x) length(x) == 2 && all(x > 0), logical(1))))
  stopifnot(length(mu_noise) == kmat, all(mu_noise >= 0))
  stopifnot(length(theta) == kmat, all(theta > 0))
  stopifnot(length(pi) == kmat, all(pi > 0 & pi < 1))
  
  if (!is.na(depth_range[1])) {
    stopifnot(is.numeric(depth_range), length(depth_range) == 2, all(depth_range > 0), depth_range[1] < depth_range[2])
  }
  
  gene_names <- rownames(B)
  cellID <- paste0("cell_", seq_len(n))
  B <- (B > 0) * 1
  A_info <- .create_adjacency_expansion(B)
  A <- A_info$A
  edges <- A_info$edge_indices
  
  matrices <- vector("list", kmat)
  for (k in seq_len(kmat)) {
    mu <- runif(p, mu_range[[k]][1], mu_range[[k]][2])
    
    edge_vals <- sample(seq(min(unlist(mu_range)), max(unlist(mu_range))), nrow(edges), replace = TRUE)
    B_weighted <- matrix(0, nrow = p, ncol = p)
    B_weighted[edges] <- edge_vals
    B_weighted <- (B_weighted | t(B_weighted)) * 1
    sigma <- B_weighted[lower.tri(B_weighted)]
    
    Y_mu <- c(mu, sigma)
    Y <- .simulate_counts_ZINB(n, Y_mu, theta[k], pi[k])
    X <- A %*% Y
    
    X <- X + .add_technical_noise(n, p, mu_noise[k], pi[k])
    X <- t(X)
    
    if (!is.null(gene_names)) colnames(X) <- gene_names
    rownames(X) <- cellID
    
    if (!is.na(depth_range[1])) {
      X <- .normalize_library_size(X, depth_range)
    }
    
    matrices[[k]] <- X
  }
  
  return(matrices)
}

