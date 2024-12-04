infer_networks <- function(count_matrices_list, method = "GENIE3", adjm = NULL, nCores = (detectCores()-1)) {
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
    
    # Loop over each matrix in the list for GENIE3, GRNBoost2, and ZILGM
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
        # Compute lambda_max and lambda_min for the current count matrix
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
      cat("Assigned lamb to lambda_results[[", z, "]] with", length(lamb), "matrices.\n")
    }
    
    return(list(network_results = network_results, lambda_results = lambda_results))
  } else {
    return(network_results)
  }
}
