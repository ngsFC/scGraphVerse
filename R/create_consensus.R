#' Create a Consensus Adjacency Matrix from Multiple Networks
#'
#' This function builds a consensus adjacency matrix from a list of binary or weighted 
#' adjacency matrices using one of three methods: \code{"vote"}, \code{"union"}, or \code{"INet"}.
#'
#' @param adj_matrix_list A list of binary adjacency matrices (square, 0/1) with identical dimensions and row/column names.
#' @param method Character string indicating the consensus strategy. One of:
#'   \itemize{
#'     \item{\code{"vote"} (default):} An edge is included if present in at least \code{threshold} proportion of matrices.
#'     \item{\code{"union"}:} An edge is included if present in at least one matrix.
#'     \item{\code{"INet"}:} Uses normalized weighted adjacency matrices and calls \code{\link[INetTool]{consensusNet}}.
#'   }
#' @param weighted_list A list of weighted adjacency matrices (required if \code{method = "INet"}).
#' @param theta Numeric. Tuning parameter passed to \code{consensusNet} (default: 0.04).
#' @param threshold Numeric in [0, 1]. Used for "vote" (fraction of matrices required to support an edge)
#'   and "INet" (thresholding inside \code{consensusNet}). Default is 0.5.
#' @param ncores Integer. Number of CPU cores to use when method is \code{"INet"}. Default is 1.
#'
#' @return A square consensus adjacency matrix (binary or weighted, depending on method).
#'
#' @details
#' The consensus is constructed differently depending on the selected method:
#' \describe{
#'   \item{\strong{vote}}{Counts edge frequency across all matrices and includes edges present in at least \code{threshold × N} matrices.}
#'   \item{\strong{union}}{Includes any edge present in any matrix.}
#'   \item{\strong{INet}}{Multiplies each binary matrix by its corresponding weighted matrix, normalizes values, and 
#'     combines networks using \code{\link[INetTool]{consensusNet}}.}
#' }
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom INetTool consensusNet
#' @export

create_consensus <- function(adj_matrix_list, method = "vote", weighted_list = NULL, theta = 0.04, threshold = 0.5, ncores = 1) {
  
  # Sum all adjacency matrices to get consensus values (default for 'vote' and 'union')
  consensus_matrix <- Reduce("+", adj_matrix_list)
  
  if (method == "vote") {
    threshold_value <- round(length(adj_matrix_list) * 0.75)
    
    # Vote threshold: 1 if >= threshold, 0 if < threshold
    consensus_matrix[consensus_matrix < threshold_value] <- 0
    consensus_matrix[consensus_matrix >= threshold_value] <- 1
  } else if (method == "union") {
    # Union method: any position with a value >= 1 in any matrix is set to 1
    consensus_matrix[consensus_matrix >= 1] <- 1
  } else if (method == "INet") {
    
    if (is.null(weighted_list)) {
      stop("For method 'INet', a weighted_list must be provided.")
    }
    
    list_weighted <- mapply(function(binary_mat, weighted_mat) {
      if (!is.matrix(binary_mat) || !is.matrix(weighted_mat)) {
        stop("Both binary and weighted elements must be matrices.")
      }
      if (!all(dim(binary_mat) == dim(weighted_mat))) {
        stop("Dimension mismatch between binary adjacency matrix and weighted matrix.")
      }
      newmat <- binary_mat * weighted_mat
      return(newmat)
    }, adj_matrix_list, weighted_list, SIMPLIFY = FALSE)
    
    # Normalize the weighted matrices
    list_genie3_norm <- lapply(list_weighted, function(mat) mat / max(mat, na.rm = TRUE))
    
    Con <- consensusNet(list_genie3_norm, theta = theta, ncores = ncores, threshold = threshold)
    consensus_matrix <- as.matrix(as_adjacency_matrix(Con$graphConsensus))
  
  } else {
    stop("Invalid method. Choose 'vote', 'union', or 'INet'.")
  }
  
  return(consensus_matrix)
}

