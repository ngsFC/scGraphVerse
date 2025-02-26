#' Infer Gene Regulatory Networks Using Different Methods
#'
#' This function infers gene regulatory networks (GRNs) from a list of count matrices using 
#' several different methods such as GENIE3, GRNBoost2, ZILGM, JRF, or PCzinb.
#'
#' @param count_matrices_list A list of matrices containing gene expression data (rows as genes and columns as samples).
#' @param method A character string specifying the method to use for inferring the network. Options are:
#'   "GENIE3", "GRNBoost2", "ZILGM", "JRF", and "PCzinb". Default is "GENIE3".
#' @param adjm An optional adjacency matrix to be used in certain methods (e.g., for PCzinb). Default is `NULL`.
#' @param nCores An integer specifying the number of cores to use for parallel computation. Default is the number of cores minus one.
#' @return A list containing the inferred networks. The exact structure depends on the method chosen.
#' @details
#' - "GENIE3" uses the GENIE3 algorithm to infer the GRN.
#' - "GRNBoost2" uses the GRNBoost2 algorithm for network inference.
#' - "ZILGM" uses the ZILGM method for network inference and also returns lambda values.
#' - "JRF" uses the JRF method, normalizing the input matrices before applying the model.
#' - "PCzinb" uses the PCzinb method for GRN inference.
#' 
#' @examples
#' \dontrun{
#' result <- infer_networks(count_matrices_list = list(matrix1, matrix2), method = "GENIE3")
#' }
#' @export
infer_networks <- function(count_matrices_list, method = "GENIE3", adjm = NULL, nCores = (BiocParallel::bpworkers(BiocParallel::bpparam()) - 1)) {
  
  # Detect and extract expression data from Seurat objects
  if (all(sapply(count_matrices_list, function(x) inherits(x, "Seurat")))) {
    message("Detected Seurat objects. Extracting expression matrices...")
    count_matrices_list <- lapply(count_matrices_list, function(obj) {
      expr_mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
      if (!inherits(expr_mat, "matrix")) {
        expr_mat <- as.matrix(expr_mat)  # Convert sparse or other formats to matrix
      }
      return(expr_mat)
    })
  }
  
  # Detect and extract logcounts from SingleCellExperiment objects
  if (all(sapply(count_matrices_list, function(x) inherits(x, "SingleCellExperiment")))) {
    message("Detected SingleCellExperiment objects. Extracting logcounts matrices...")
    count_matrices_list <- lapply(count_matrices_list, function(sce) {
      expr_mat <- SummarizedExperiment::assay(sce, "logcounts")
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
  
  network_results <- list()
  
  # Apply JRF with normalization (Parallelized with BiocParallel)
  if (method == "JRF") {
    count_matrices_list <- BiocParallel::bplapply(count_matrices_list, function(x) {
      (x - mean(x)) / sd(x)
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))
    
    netout <- JRF::JRF(X = count_matrices_list, 
                       genes.name = rownames(count_matrices_list[[1]]), 
                       ntree = 500, 
                       mtry = round(sqrt(length(rownames(count_matrices_list[[1]])) - 1)))
    
    network_results[[1]] <- netout
    return(network_results)
  }
  
  # Apply methods that handle parallelization internally
  if (method %in% c("GENIE3", "ZILGM", "PCzinb")) {
    network_results <- lapply(count_matrices_list, function(count_matrix) {
      if (method == "GENIE3") {
        return(GENIE3::getLinkList(GENIE3::GENIE3(count_matrix, nCores = nCores)))
      } else if (method == "ZILGM") {
        lambda_max <- ZILGM::find_lammax(t(as.matrix(count_matrix)))
        lambda_min <- 1e-4 * lambda_max
        lambs <- exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
        nb2_fit <- ZILGM::zilgm(X = t(as.matrix(count_matrix)), lambda = lambs, family = "NBII", update_type = "IRLS", do_boot = TRUE, boot_num = 10, sym = "OR", nCores = nCores)
        return(nb2_fit$network[[nb2_fit$opt_index]])
      } else if (method == "PCzinb") {
        netout <- learn2count::PCzinb(t(as.matrix(count_matrix)), method = "zinb1", maxcard = 2)
        rownames(netout) <- rownames(adjm)
        colnames(netout) <- colnames(adjm)
        return(netout)
      }
    })
    return(network_results)
  }
  
  # Apply BiocParallel only for GRNBoost2
  if (method == "GRNBoost2") {
    network_results <- BiocParallel::bplapply(seq_along(count_matrices_list), function(j) {
      count_matrix <- count_matrices_list[[j]]
      count_matrix_df <- as.data.frame(t(count_matrix))
      genes <- colnames(count_matrix_df)
      df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), columns = genes, index = rownames(count_matrix_df))
      return(arboreto$grnboost2(df_pandas, gene_names = genes))
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  }
  
  return(network_results)
}
