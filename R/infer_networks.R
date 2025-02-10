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
  library(BiocParallel)
  
  # Detect and extract expression data from Seurat objects
  if (all(sapply(count_matrices_list, function(x) inherits(x, "Seurat")))) {
    message("Detected Seurat objects. Extracting expression matrices...")
    count_matrices_list <- bplapply(count_matrices_list, function(obj) {
      as.matrix(t(GetAssayData(obj, assay = "RNA", layer = "data")))
    }, BPPARAM = MulticoreParam(nCores))
  }
  
  # Detect and extract logcounts from SingleCellExperiment objects
  if (all(sapply(count_matrices_list, function(x) inherits(x, "SingleCellExperiment")))) {
    message("Detected SingleCellExperiment objects. Extracting logcounts matrices...")
    count_matrices_list <- bplapply(count_matrices_list, function(sce) {
      as.matrix(t(assay(sce, "logcounts")))
    }, BPPARAM = MulticoreParam(nCores))
  }
  
  # Validate method input
  valid_methods <- c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  
  network_results <- list()
  
  # Apply JRF with normalization
  if (method == "JRF") {
    count_matrices_list <- bplapply(count_matrices_list, function(x) {
      (x - mean(x)) / sd(x)
    }, BPPARAM = MulticoreParam(nCores))
    
    netout <- JRF(X = count_matrices_list, 
                  genes.name = rownames(count_matrices_list[[1]]), 
                  ntree = 500, 
                  mtry = round(sqrt(length(rownames(count_matrices_list[[1]])) - 1)))
    
    network_results[[1]] <- netout
    return(network_results)
  }
  
  # Parallel network inference for other methods
  network_results <- bplapply(seq_along(count_matrices_list), function(j) {
    count_matrix <- count_matrices_list[[j]]
    
    if (method == "GENIE3") {
      return(getLinkList(GENIE3(t(count_matrix), nCores = nCores)))
    } else if (method == "GRNBoost2") {
      count_matrix_df <- as.data.frame(count_matrix)
      genes <- colnames(count_matrix_df)
      df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), columns = genes, index = rownames(count_matrix_df))
      return(arboreto$grnboost2(df_pandas, gene_names = genes))
    } else if (method == "ZILGM") {
      lambda_max <- find_lammax(as.matrix(count_matrix))
      lambda_min <- 1e-4 * lambda_max
      lambs <- exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
      nb2_fit <- zilgm(X = as.matrix(count_matrix), lambda = lambs, family = "NBII", update_type = "IRLS", do_boot = TRUE, boot_num = 10, sym = "OR", nCores = nCores)
      return(nb2_fit$network[[nb2_fit$opt_index]])
    } else if (method == "PCzinb") {
      netout <- PCzinb(as.matrix(count_matrix), method = "zinb1", maxcard = 2)
      rownames(netout) <- rownames(adjm)
      colnames(netout) <- colnames(adjm)
      return(netout)
    }
  }, BPPARAM = MulticoreParam(nCores))
  
  return(network_results)
}
