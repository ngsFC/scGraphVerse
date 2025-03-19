# Helper function to convert a link list into a weighted adjacency matrix.
# Assumes link_list has columns: regulatoryGene, targetGene, weight.
links_to_adjacency <- function(link_list, genes) {
  n <- length(genes)
  adj <- matrix(0, n, n)
  rownames(adj) <- genes
  colnames(adj) <- genes
  for(i in seq_len(nrow(link_list))) {
    reg <- as.character(link_list[i, "regulatoryGene"])
    tgt <- as.character(link_list[i, "targetGene"])
    wt  <- link_list[i, "weight"]
    if(reg %in% genes && tgt %in% genes) {
      adj[reg, tgt] <- wt
    }
  }
  return(adj)
}

# Provided symmetrize function.
symmetrize <- function(matrix_list, weight_function = "mean", nCores = BiocParallel::bpworkers(BiocParallel::bpparam())) {
  if (!is.list(matrix_list) || !all(sapply(matrix_list, is.matrix))) {
    stop("matrix_list must be a list of matrices")
  }
  
  weight_function <- match.fun(weight_function)  # Ensure valid function input
  
  # Parallel processing of matrix symmetrization
  symmetrized_matrices <- BiocParallel::bplapply(matrix_list, function(mat) {
    p <- nrow(mat)
    sym_mat <- mat  # Copy structure
    
    for (i in seq_len(p - 1)) {
      for (j in seq(i + 1, p)) {
        val_ij <- mat[i, j]
        val_ji <- mat[j, i]
        
        if (val_ij == 0 || val_ji == 0) {
          symmetric_value <- max(val_ij, val_ji)  # Use the non-zero value
        } else {
          symmetric_value <- weight_function(c(val_ij, val_ji))  # Apply function
        }
        
        sym_mat[i, j] <- symmetric_value
        sym_mat[j, i] <- symmetric_value
      }
    }
    return(sym_mat)
  }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  
  return(symmetrized_matrices)
}

#' Infer Gene Regulatory Networks with Stability Selection and Symmetrization
#'
#' This function infers gene regulatory networks (GRNs) from a list of count matrices using 
#' several different methods (e.g., GENIE3, GRNBoost2, ZILGM, JRF, or PCzinb) with an integrated 
#' stability selection procedure. For each count matrix, many subsamples are drawn, the network is 
#' inferred on each subsample and thresholded into a binary adjacency matrix (using a 99th percentile 
#' cutoff on the weighted network). Then, the selection frequency is computed over all subsamples, 
#' and a final binary, symmetric adjacency matrix is returned that retains only those edges that appear 
#' with a frequency above a specified threshold.
#'
#' @param count_matrices_list A list of matrices (or Seurat/SingleCellExperiment objects) containing gene expression data (rows as genes and columns as samples).
#' @param method A character string specifying the method to use for inferring the network. Options are:
#'   "GENIE3", "GRNBoost2", "ZILGM", "JRF", and "PCzinb". Default is "GENIE3".
#' @param adjm An optional adjacency matrix (used for some methods such as PCzinb). Default is NULL.
#' @param nCores An integer specifying the number of cores to use for parallel computation. Default is the number of cores minus one.
#' @param n_subsamples Number of subsamples (bootstrap iterations) for stability selection. Default is 50.
#' @param subsample_fraction Fraction of the samples (columns) to use for each subsample. Default is 0.7.
#' @param threshold_frequency Selection frequency threshold (between 0 and 1) for an edge to be kept. Default is 0.7.
#' @return A list of binary, symmetric adjacency matrices (one per count matrix) after applying stability selection.
#' @export
infer_networks <- function(count_matrices_list, 
                           method = "GENIE3", 
                           adjm = NULL, 
                           nCores = (BiocParallel::bpworkers(BiocParallel::bpparam()) - 1),
                           n_subsamples = 50,
                           subsample_fraction = 0.7,
                           threshold_frequency = 0.7) {
  
  # --- Preprocess count matrices (Seurat / SingleCellExperiment) ---
  if (all(sapply(count_matrices_list, function(x) inherits(x, "Seurat")))) {
    message("Detected Seurat objects. Extracting expression matrices...")
    count_matrices_list <- lapply(count_matrices_list, function(obj) {
      expr_mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
      if (!inherits(expr_mat, "matrix")) {
        expr_mat <- as.matrix(expr_mat)
      }
      return(expr_mat)
    })
  }
  
  if (all(sapply(count_matrices_list, function(x) inherits(x, "SingleCellExperiment")))) {
    message("Detected SingleCellExperiment objects. Extracting counts matrices...")
    count_matrices_list <- lapply(count_matrices_list, function(sce) {
      expr_mat <- SummarizedExperiment::assay(sce, "counts")
      if (!inherits(expr_mat, "matrix")) {
        expr_mat <- as.matrix(expr_mat)
      }
      return(expr_mat)
    })
  }
  
  # Validate method input
  valid_methods <- c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  
  # --- Helper function: infer network on a subsample ---
  infer_network_single <- function(count_matrix_subsample) {
    if (method == "GENIE3") {
      links <- GENIE3::GENIE3(count_matrix_subsample, nCores = nCores)
      link_list <- GENIE3::getLinkList(links)
      genes <- rownames(count_matrix_subsample)
      weighted_adj <- links_to_adjacency(link_list, genes)
      weighted_adj <- (weighted_adj + t(weighted_adj)) / 2  # Symmetrize weighted matrix
      return(weighted_adj)
      
    } else if (method == "GRNBoost2") {
      count_matrix_df <- as.data.frame(t(count_matrix_subsample))
      genes <- colnames(count_matrix_df)
      df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), 
                                    columns = genes, 
                                    index = rownames(count_matrix_df))
      result <- arboreto$grnboost2(df_pandas, gene_names = genes)
      link_list <- as.data.frame(result)
      weighted_adj <- links_to_adjacency(link_list, rownames(count_matrix_subsample))
      weighted_adj <- (weighted_adj + t(weighted_adj)) / 2
      return(weighted_adj)
      
    } else if (method == "ZILGM") {
      lambda_max <- ZILGM::find_lammax(t(as.matrix(count_matrix_subsample)))
      lambda_min <- 1e-4 * lambda_max
      lambs <- exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
      nb2_fit <- ZILGM::zilgm(X = t(as.matrix(count_matrix_subsample)), 
                              lambda = lambs, 
                              family = "NBII", 
                              update_type = "IRLS", 
                              do_boot = TRUE, 
                              boot_num = 10, 
                              sym = "OR", 
                              nCores = nCores)
      return(nb2_fit$network[[nb2_fit$opt_index]])
      
    } else if (method == "PCzinb") {
      netout <- learn2count::PCzinb(t(as.matrix(count_matrix_subsample)), method = "zinb1", maxcard = 2)
      if (!is.null(adjm)) {
        rownames(netout) <- rownames(adjm)
        colnames(netout) <- colnames(adjm)
      }
      return(netout)
      
    } else if (method == "JRF") {
      norm_matrix <- (count_matrix_subsample - mean(count_matrix_subsample)) / sd(count_matrix_subsample)
      netout <- JRF::JRF(X = list(norm_matrix), 
                         genes.name = rownames(norm_matrix), 
                         ntree = 500, 
                         mtry = round(sqrt(length(rownames(norm_matrix)) - 1)))
      return(netout)
    }
  }
  
  # --- Helper function: threshold weighted adjacency into binary ---
  threshold_weighted_adj <- function(weighted_adj) {
    # Compute the 99th percentile from the upper triangular part of the matrix.
    weights <- weighted_adj[upper.tri(weighted_adj)]
    cutoff <- quantile(weights, 0.99, na.rm = TRUE)
    binary_adj <- ifelse(weighted_adj > cutoff, 1, 0)
    return(binary_adj)
  }
  
  # --- Stability Selection: Loop over each count matrix ---
  stability_selected_networks <- lapply(count_matrices_list, function(count_matrix) {
    genes <- rownames(count_matrix)
    selection_counts <- matrix(0, nrow = length(genes), ncol = length(genes))
    rownames(selection_counts) <- genes
    colnames(selection_counts) <- genes
    
    subsample_results <- BiocParallel::bplapply(seq_len(n_subsamples), function(i) {
      subsample_indices <- sample(seq_len(ncol(count_matrix)), 
                                  size = round(ncol(count_matrix) * subsample_fraction))
      subsample <- count_matrix[, subsample_indices, drop = FALSE]
      
      weighted_adj <- infer_network_single(subsample)
      binary_adj <- threshold_weighted_adj(weighted_adj)
      return(binary_adj)
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))
    
    # Sum binary selections over all subsamples.
    for (binary_adj in subsample_results) {
      selection_counts <- selection_counts + binary_adj
    }
    
    # Compute selection frequency for each edge.
    selection_frequency <- selection_counts / n_subsamples
    
    # Final binary network: retain edges with frequency at least threshold_frequency.
    final_binary_network <- ifelse(selection_frequency >= threshold_frequency, 1, 0)
    
    # Ensure the final network is symmetric.
    final_binary_network <- symmetrize(list(final_binary_network), weight_function = "mean", nCores = nCores)[[1]]
    return(final_binary_network)
  })
  
  return(stability_selected_networks)
}

