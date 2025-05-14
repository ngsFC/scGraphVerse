#' Simulate Zero-Inflated Negative Binomial (ZINB) Count Matrices with
#' Sequencing Depth
#'
#' Simulates one or more count matrices following a zero-inflated negative
#' binomial (ZINB) distribution, incorporating gene-gene interaction
#' structures and cell-specific sequencing depth variation.
#'
#' @param n Integer. Number of cells (samples) in each simulated matrix.
#' @param p Integer. Number of genes (features) in each simulated matrix.
#' @param B A symmetric binary adjacency matrix (0/1) defining gene-gene
#'   connectivity. Row and column names correspond to gene names.
#' @param mu_range List of numeric vectors (length 2 each). Range of gene
#'   expression means for each simulated matrix.
#' @param mu_noise Numeric vector. Mean of background noise for each matrix.
#' @param theta Numeric vector. Dispersion parameters of the negative
#'   binomial distribution for each matrix. Smaller \code{theta} implies
#'   higher overdispersion.
#' @param pi Numeric vector. Probability of excess zeros (\code{0 < pi < 1})
#'   for each matrix.
#' @param kmat Integer. Number of count matrices to simulate. Default is 1.
#' @param depth_range Numeric vector of length 2 or \code{NA}. Range of total
#'   sequencing depth per cell. If \code{NA}, no depth adjustment is
#'   performed.
#'
#' @return A list containing \code{kmat} matrices. Each matrix has:
#'   \itemize{
#'     \item Rows representing cells (\code{cell_1, ..., cell_n}).
#'     \item Columns representing genes (\code{rownames(B)}).
#'     \item Count values following a ZINB distribution.
#'   }
#'
#' @details Each simulated matrix:
#'   \enumerate{
#'     \item Generates gene expression values based on a ZINB model.
#'     \item Modulates expression using the adjacency matrix \code{B}.
#'     \item Applies random sequencing depth scaling if
#'       \code{depth_range} is provided.
#'   }
#'
#' Useful for benchmarking single-cell RNA-seq network inference methods
#' with dropout events and network structure.
#'
#' @importFrom distributions3 rzinbinom
#' @export
#'
#' @examples
#' data(adj_truth)
#' nodes <- nrow(adj_truth)
#' sims <- zinb_simdata(
#'     n = 50,
#'     p = nodes,
#'     B = adj_truth,
#'     mu_range = list(c(1, 4), c(1, 7), c(1, 10)),
#'     mu_noise = c(1, 3, 5),
#'     theta = c(1, 0.7, 0.5),
#'     pi = c(0.2, 0.2, 0.2),
#'     kmat = 3,
#'     depth_range = c(0.8 * nodes * 3, 1.2 * nodes * 3)
#' )
zinb_simdata <- function(
    n,
    p,
    B,
    mu_range,
    mu_noise,
    theta,
    pi,
    kmat = 1,
    depth_range = NA) {
    stopifnot(
        is.numeric(n), n > 0, floor(n) == n,
        is.numeric(p), p > 0, floor(p) == p,
        is.matrix(B), nrow(B) == ncol(B), all(B %in% c(0, 1)),
        is.numeric(kmat), kmat > 0, floor(kmat) == kmat,
        length(mu_range) == kmat,
        all(vapply(mu_range, function(x) {
            length(x) == 2 && all(x > 0)
        }, logical(1))),
        length(mu_noise) == kmat, all(mu_noise >= 0),
        length(theta) == kmat, all(theta > 0),
        length(pi) == kmat, all(pi > 0 & pi < 1)
    )

    if (!is.na(depth_range[1])) {
        stopifnot(
            is.numeric(depth_range),
            length(depth_range) == 2,
            all(depth_range > 0),
            depth_range[1] < depth_range[2]
        )
    }

    gene_names <- rownames(B)
    cellID <- paste0("cell_", seq_len(n))
    B <- (B > 0) * 1
    A_info <- .create_adjacency_expansion(B)
    A <- A_info$A
    edges <- A_info$edge_indices

    B[edges] <- sample(
        seq(min(unlist(mu_range)), max(unlist(mu_range))),
        length(edges[, 1]),
        replace = TRUE
    )
    B <- (B | t(B)) * 1

    matrices <- vector("list", kmat)
    for (k in seq_len(kmat)) {
        mu <- runif(p, mu_range[[k]][1], mu_range[[k]][2])
        sigma <- B
        nonzero_sigma <- sigma[lower.tri(sigma) & sigma != 0]
        Y_mu <- c(mu, nonzero_sigma)
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

    matrices
}
