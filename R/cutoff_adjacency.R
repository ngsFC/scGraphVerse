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
                             total_cores = BiocParallel::bpworkers(BiocParallel::bpparam()),
                             grnboost_modules = NULL,
                             seed = 123,
                             debug = FALSE) {
  
  method <- match.arg(method, choices = c("GENIE3", "GRNBoost2", "JRF"))
  weight_function <- match.fun(weight_function)
  
  if (length(count_matrices) != length(weighted_adjm_list)) {
    stop("Length of count_matrices must match weighted_adjm_list.")
  }
  
  # Prepare expression matrices
  count_matrices <- lapply(count_matrices, function(obj) {
    if (inherits(obj, "Seurat")) {
      as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
    } else if (inherits(obj, "SingleCellExperiment")) {
      as.matrix(SummarizedExperiment::assay(obj, "counts"))
    } else {
      as.matrix(obj)
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
  
  create_shuffled_matrices <- function(mat, n, base_seed) {
    lapply(seq_len(n), function(i) {
      seed_vector <- sample(base_seed + i + seq_len(nrow(mat)))
      shuffle_rows(mat, seed_vector)
    })
  }
  
  # Define total jobs (k matrices × n shuffles)
  job_list <- expand.grid(matrix_idx = seq_along(count_matrices), shuffle_idx = seq_len(n))
  job_list <- lapply(seq_len(nrow(job_list)), function(i) {
    list(matrix_idx = job_list$matrix_idx[i], shuffle_idx = job_list$shuffle_idx[i])
  })
  
  # Parallel strategy
  nCores_outer <- min(total_cores, length(job_list))
  nCores_inner <- max(floor(total_cores / nCores_outer), 1)
  
  param_outer <- if (method == "GRNBoost2") {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::MulticoreParam(workers = nCores_outer, RNGseed = seed)
  }
  
  percentile_values_by_matrix <- vector("list", length(count_matrices))
  
  results <- BiocParallel::bplapply(job_list, function(job) {
    mat_idx <- job$matrix_idx
    shuf_idx <- job$shuffle_idx
    mat <- count_matrices[[mat_idx]]
    
    base_seed <- as.integer(round(seed * 1000 + mat_idx * 10 + shuf_idx))
    shuffled <- create_shuffled_matrices(mat, 1, base_seed)[[1]]
    
    args <- list(count_matrices_list = list(shuffled), method = method, adjm = NULL, seed = base_seed)
    
    if (method == "GENIE3") {
      args$total_cores <- nCores_inner
    }
    if (method == "GRNBoost2") {
      args$grnboost_modules <- grnboost_modules
      args$total_cores <- nCores_inner
    }
    
    inferred <- do.call(infer_networks, args)
    adjm <- generate_adjacency(inferred)
    symm <- symmetrize(adjm, weight_function = weight_function)[[1]]
    q_value <- quantile(symm[upper.tri(symm)], quantile_threshold, names = FALSE)
    list(matrix_idx = mat_idx, q_value = q_value)
  }, BPPARAM = param_outer)
  
  # Aggregate results by matrix
  for (res in results) {
    mat_idx <- res$matrix_idx
    percentile_values_by_matrix[[mat_idx]] <- c(percentile_values_by_matrix[[mat_idx]], res$q_value)
  }
  
  binary_adjm_list <- lapply(seq_along(weighted_adjm_list), function(idx) {
    avg_cutoff <- mean(percentile_values_by_matrix[[idx]])
    if (debug) message(sprintf("[Method: %s] Matrix %d → Cutoff = %.5f", method, idx, avg_cutoff))
    ifelse(weighted_adjm_list[[idx]] > avg_cutoff, 1, 0)
  })
  
  return(binary_adjm_list)
}
