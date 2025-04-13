#' Apply Cutoff to Adjacency Matrices Based on Shuffled Network Percentiles
#'
#' Thresholds weighted adjacency matrices by estimating cutoffs from shuffled networks.
#' Supports: "GENIE3", "GRNBoost2", "JRF".
#'
#' @param count_matrices A list of expression matrices (genes × cells) or Seurat/SCE objects.
#' @param weighted_adjm_list A list of weighted adjacency matrices to threshold.
#' @param n Integer. Number of shuffled matrices to generate per original matrix.
#' @param method Character. One of: "GENIE3", "GRNBoost2", "JRF".
#' @param quantile_threshold Numeric. Default is 0.99. Quantile used for cutoff estimation.
#' @param weight_function Character or function for symmetrizing adjacency matrices (e.g., "mean", "max").
#' @param nCores Number of CPU cores for GENIE3 or JRF parallelism.
#' @param grnboost_modules Python modules for GRNBoost2.
#' @param seed Optional seed for reproducibility.
#' @param debug Logical. If TRUE, prints progress and cutoff info.
#'
#' @return A list of binary adjacency matrices (same dimension as input), thresholded by estimated cutoffs.
#' @export
cutoff_adjacency <- function(count_matrices, 
                             weighted_adjm_list, 
                             n, 
                             method = "GENIE3", 
                             quantile_threshold = 0.99,
                             weight_function = "mean", 
                             nCores = 1,
                             grnboost_modules = NULL,
                             seed = NULL,
                             debug = FALSE) {
  
  method <- match.arg(method, choices = c("GENIE3", "GRNBoost2", "JRF"))
  weight_function <- match.fun(weight_function)
  
  if (length(count_matrices) != length(weighted_adjm_list)) {
    stop("Length of count_matrices must match weighted_adjm_list.")
  }
  
  count_matrices <- lapply(count_matrices, function(obj) {
    if (inherits(obj, "Seurat")) {
      return(as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")))
    } else if (inherits(obj, "SingleCellExperiment")) {
      return(as.matrix(SummarizedExperiment::assay(obj, "counts")))
    } else {
      return(as.matrix(obj))
    }
  })
  
  shuffle_rows <- function(matrix, seed_vector) {
    shuffled <- matrix
    for (i in seq_len(nrow(matrix))) {
      set.seed(seed_vector[i])
      shuffled[i, ] <- sample(matrix[i, ])
    }
    shuffled
  }
  
  create_shuffled_matrices <- function(mat, n, base_seed = 1234) {
    lapply(seq_len(n), function(i) {
      seed_vector <- sample(base_seed + i + seq_len(nrow(mat)))
      shuffle_rows(mat, seed_vector)
    })
  }
  
  binary_adjm_list <- lapply(seq_along(count_matrices), function(idx) {
    mat <- count_matrices[[idx]]
    shuffled_list <- create_shuffled_matrices(mat, n, base_seed = if (!is.null(seed)) seed else 1000 + idx)
    
    percentile_values <- lapply(seq_along(shuffled_list), function(shuf_idx) {
      shuf_mat <- shuffled_list[[shuf_idx]]
      args <- list(count_matrices_list = list(shuf_mat), method = method, adjm = NULL, seed = seed)
      
      if (method == "GENIE3") {
        args$nCores <- nCores
      }
      
      if (method == "GRNBoost2") {
        args$grnboost_modules <- grnboost_modules
      }
      
      if (method == "JRF") {
        args$nCores <- nCores
        inferred <- BiocParallel::bplapply(
          list(shuf_mat),
          function(x) do.call(infer_networks, args),
          BPPARAM = BiocParallel::MulticoreParam(nCores)
        )[[1]]
      } else {
        inferred <- do.call(infer_networks, args)
      }
      
      adjm <- generate_adjacency(inferred)
      symm <- symmetrize(adjm, weight_function = weight_function)[[1]]
      quantile(symm[upper.tri(symm)], quantile_threshold, names = FALSE)
      
      # Adding delay and garbage collection to handle GRNBoost2 Dask port conflicts
      if (method == "GRNBoost2") {
        Sys.sleep(5)
        gc()
      }
    })
    
    avg_cutoff <- mean(unlist(percentile_values))
    if (debug) message(sprintf("[Method: %s] Matrix %d → Cutoff = %.5f", method, idx, avg_cutoff))
    
    binary_mat <- ifelse(weighted_adjm_list[[idx]] > avg_cutoff, 1, 0)
    binary_mat
  })
  
  return(binary_adjm_list)
}
