#' Infer Gene Regulatory Networks Using Multiple Methods
#'
#' This function infers gene regulatory networks (GRNs) from a list of gene expression matrices 
#' using various methods such as GENIE3, GRNBoost2, ZILGM, JRF, or PCzinb.
#'
#' @param count_matrices_list A list of matrices containing gene expression data (rows as genes, columns as samples).
#'   Can also contain Seurat or SingleCellExperiment objects, which will be converted to count matrices.
#' @param method A character string specifying the method for network inference. Options include:
#'   \itemize{
#'     \item "GENIE3" - Uses the GENIE3 algorithm to infer the GRN.
#'     \item "GRNBoost2" - Uses the GRNBoost2 algorithm.
#'     \item "ZILGM" - Uses the ZILGM model, returning a sparse adjacency matrix.
#'     \item "JRF" - Uses the JRF method with normalized inputs.
#'     \item "PCzinb" - Uses the PCzinb method (requires an adjacency matrix `adjm`).
#'   }
#'   Default is "GENIE3".
#' @param adjm An optional adjacency matrix required for the PCzinb method. Default is `NULL`.
#' @param nCores An integer specifying the number of cores to use for parallel computation. 
#'   Default is the number of available cores minus one.
#'
#' @return A list containing inferred regulatory networks. The structure of the output depends on the method:
#'   \itemize{
#'     \item "GENIE3" - A data frame of regulatory links with columns: `regulatoryGene`, `targetGene`, and `weight`.
#'     \item "GRNBoost2" - A list of inferred networks per dataset.
#'     \item "ZILGM" - A sparse adjacency matrix representing inferred connections.
#'     \item "JRF" - A fitted JRF model object.
#'     \item "PCzinb" - A sparse inferred network.
#'   }
#'
#' @examples
#' \dontrun{
#' result <- infer_networks(count_matrices_list = list(matrix1, matrix2), method = "GENIE3")
#' }
#' @export
infer_networks <- function(count_matrices_list, method = "GENIE3", adjm = NULL, 
                           nCores = (BiocParallel::bpworkers(BiocParallel::bpparam()) - 1)) {
  
  # Validate method input
  valid_methods <- c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb")
  method <- match.arg(method, choices = valid_methods)
  
  # Convert different input formats to numeric matrices
  count_matrices_list <- lapply(count_matrices_list, function(obj) {
    if (inherits(obj, "Seurat")) {
      expr_mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
    } else if (inherits(obj, "SingleCellExperiment")) {
      expr_mat <- SummarizedExperiment::assay(obj, "counts")
    } else if (is.matrix(obj) || is.data.frame(obj)) {
      expr_mat <- as.matrix(obj)
    } else {
      stop("Invalid input type. Provide matrices, Seurat, or SingleCellExperiment objects.")
    }
    
    # Ensure matrix format
    if (!inherits(expr_mat, "matrix")) {
      expr_mat <- Matrix::as.matrix(expr_mat)  
    }
    
    # Ensure numeric data
    if (!is.numeric(expr_mat)) {
      stop("Count matrices must contain numeric values.")
    }
    
    return(expr_mat)
  })
  
  network_results <- list()

  # GENIE3: Return edge list instead of adjacency matrix
  if (method == "GENIE3") {
    network_results <- lapply(count_matrices_list, function(count_matrix) {
      adj_matrix <- GENIE3::GENIE3(count_matrix, nCores = nCores)
      return(GENIE3::getLinkList(adj_matrix))  # Converts to edge list
    })
    return(network_results)
  }
  
  # GRNBoost2: Ensure Python conversion works properly
  if (method == "GRNBoost2") {
    network_results <- BiocParallel::bplapply(seq_along(count_matrices_list), function(j) {
      count_matrix <- count_matrices_list[[j]]
      count_matrix_df <- as.data.frame(t(count_matrix))  # Convert to DataFrame
      genes <- colnames(count_matrix_df)
      
      df_pandas <- reticulate::r_to_py(count_matrix_df)  # Convert to Python Pandas
      return(arboreto$grnboost2(df_pandas, gene_names = genes))
      
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))
    return(network_results)
  }
  
  # ZILGM Method
  if (method == "ZILGM") {
    network_results <- lapply(count_matrices_list, function(count_matrix) {
      lambda_max <- ZILGM::find_lammax(t(count_matrix))
      lambda_min <- 1e-4 * lambda_max
      lambs <- exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
      nb2_fit <- ZILGM::zilgm(X = t(count_matrix), lambda = lambs, family = "NBII", 
                              update_type = "IRLS", do_boot = TRUE, boot_num = 10, 
                              sym = "OR", nCores = nCores)
      return(nb2_fit$network[[nb2_fit$opt_index]])
    })
    return(network_results)
  }

  # JRF: Apply normalization before inference
  if (method == "JRF") {
    count_matrices_list <- BiocParallel::bplapply(count_matrices_list, function(x) {
      apply(x, 1, function(g) (g - mean(g)) / sd(g))
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))
    
    netout <- JRF::JRF(X = count_matrices_list, 
                       genes.name = rownames(count_matrices_list[[1]]), 
                       ntree = 500, 
                       mtry = round(sqrt(nrow(count_matrices_list[[1]]))))
    network_results[[1]] <- netout
    return(network_results)
  }

  # PCzinb Method (Requires adjacency matrix)
  if (method == "PCzinb") {
    if (is.null(adjm)) stop("PCzinb method requires a predefined adjacency matrix (adjm).")
    
    network_results <- lapply(count_matrices_list, function(count_matrix) {
      netout <- learn2count::PCzinb(t(count_matrix), method = "zinb1", maxcard = 2)
      rownames(netout) <- rownames(adjm)
      colnames(netout) <- colnames(adjm)
      return(netout)
    })
    return(network_results)
  }

  return(network_results)
}
