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

infer_networks <- function(count_matrices_list, method = "GENIE3", adjm = NULL, nCores = (detectCores() - 1)) {
  # Detect if input is a list of Seurat objects and extract expression matrices
  if (all(sapply(count_matrices_list, function(x) inherits(x, "Seurat")))) {
    message("Detected Seurat objects. Extracting expression matrices...")
    
    count_matrices_list <- lapply(count_matrices_list, function(obj) {
      GetAssayData(obj, assay = "RNA", layer = "data") # Extract expression data
    })
    
    count_matrices_list <- lapply(count_matrices_list, as.matrix) # Convert to matrix format
    count_matrices_list <- lapply(count_matrices_list, t) # Convert to matrix format
  }
  
  # Detect if input is a list of SingleCellExperiment (SCE) objects and extract logcounts
  if (all(sapply(count_matrices_list, function(x) inherits(x, "SingleCellExperiment")))) {
    message("Detected SingleCellExperiment objects. Extracting logcounts matrices...")
    
    count_matrices_list <- lapply(count_matrices_list, function(sce) {
      assay(sce, "logcounts") # Extract log-normalized counts
    })
    
    count_matrices_list <- lapply(count_matrices_list, as.matrix) # Convert to matrix format
    count_matrices_list <- lapply(count_matrices_list, t) # Convert to matrix format
  }
  
  # Validate method input
  if (!method %in% c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb")) {
    stop("Invalid method. Choose either 'GENIE3', 'GRNBoost2', 'ZILGM', 'PCzinb', or 'JRF'.")
  }
  
  network_results <- list()
  
  if (method == "JRF") {
    # Normalize matrices for JRF
    jrf_matrices <- lapply(count_matrices_list, t)
    jrf_matrices <- lapply(jrf_matrices, function(x) {
      (x - mean(x)) / sd(x)
    })
    
    netout <- JRF(X = jrf_matrices, 
                  genes.name = rownames(jrf_matrices[[1]]), 
                  ntree = 500, 
                  mtry = round(sqrt(length(rownames(jrf_matrices[[1]])) - 1)))
    
    network_results[[1]] <- netout
    
  } else {
    mlamb <- list()
    llamb <- list()
    
    # Loop over each matrix for GENIE3, GRNBoost2, and ZILGM
    for (j in seq_along(count_matrices_list)) {
      if (method == "GENIE3") {
        regulatory_network <- GENIE3(t(count_matrices_list[[j]]), nCores = nCores)
        regulatory_network <- getLinkList(regulatory_network)
        netout <- regulatory_network
        network_results[[j]] <- netout
        
      } else if (method == "GRNBoost2") {
        count_matrix_df <- as.data.frame(count_matrices_list[[j]])
        genes <- colnames(count_matrix_df)
        
        df_pandas <- pandas$DataFrame(data = as.matrix(count_matrix_df), 
                                      columns = genes, 
                                      index = rownames(count_matrix_df))
        
        netout <- arboreto$grnboost2(df_pandas, gene_names = genes)
        network_results[[j]] <- netout
        
      } else if (method == "ZILGM") {
        # Compute lambda_max and lambda_min for ZILGM
        lambda_max <- find_lammax(as.matrix(count_matrices_list[[j]]))
        lambda_min <- 1e-4 * lambda_max
        lambs <- exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
        
        nb2_fit <- zilgm(X = as.matrix(count_matrices_list[[j]]), lambda = lambs, family = "NBII",
                         update_type = "IRLS", do_boot = TRUE, boot_num = 10, sym = "OR", nCores = nCores)
        
        netout <- nb2_fit$network[[nb2_fit$opt_index]]
        network_results[[j]] <- netout
        llamb[[j]] <- nb2_fit$lambda
        mlamb[[j]] <- nb2_fit$network
        
      } else if (method == "PCzinb") {
        netout <- PCzinb(as.matrix(count_matrices_list[[j]]), method="zinb1", maxcard=2)
        rownames(netout) <- rownames(adjm)
        colnames(netout) <- colnames(adjm)
        network_results[[j]] <- netout
      }
    }
  }
  
  # Process results for ZILGM to include lambda results
  if (method == "ZILGM") {
    lambda_results <- vector("list", length(mlamb))
    
    network_results <- lapply(network_results, as.matrix)
    if (!is.null(adjm)) {
      for (k in seq_along(network_results)) {
        rownames(network_results[[k]]) <- rownames(adjm)
        colnames(network_results[[k]]) <- colnames(adjm)
      }
    }
    
    for (z in seq_along(mlamb)) {
      lamb <- mlamb[[z]]
      lamb <- lapply(lamb, as.matrix)
      names(lamb) <- llamb[[z]]
      
      if (!is.null(adjm)) {
        for (k in seq_along(lamb)) {
          rownames(lamb[[k]]) <- rownames(adjm)
          colnames(lamb[[k]]) <- colnames(adjm)
        }
      }
      
      lambda_results[[z]] <- lamb
      message("Assigned lambda results for dataset ", z, " with ", length(lamb), " matrices.")
    }
    
    return(list(network_results = network_results, lambda_results = lambda_results))
  } else {
    return(network_results)
  }
}
