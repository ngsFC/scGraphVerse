#' Generate Adjacency Matrices from a List of Data Frames
#'
#' This function generates adjacency matrices based on the information provided in 
#' a list of data frames. Each data frame is expected to have three columns: 
#' the first and second columns represent gene pairs, and the third column represents 
#' the weight of the interaction between the genes.
#'
#' @param df_list A list of data frames, each containing three columns: 
#'   - First column: Gene1 (character)
#'   - Second column: Gene2 (character)
#'   - Third column: Weight (numeric)
#' @return A list of adjacency matrices, one for each data frame in `df_list`.
#' @details 
#' The function creates an adjacency matrix for each data frame in `df_list` by placing 
#' the interaction weights in the correct positions based on the gene pairings. The 
#' diagonal of the resulting matrices is set to zero to ensure no self-interactions are included.
#' 
#' @examples
#' 
#' # Example data frame
#' df1 <- data.frame(Gene1 = c("GeneA", "GeneB"), Gene2 = c("GeneB", "GeneC"), Weight = c(0.5, 0.8))
#' df_list <- list(df1)
#' adjacency_matrices <- generate_adjacency(df_list)
#' 
#' @export

generate_adjacency <- function(df_list, nCores = BiocParallel::bpworkers(BiocParallel::bpparam())) {
  if (!is.list(df_list) || !all(sapply(df_list, is.data.frame))) {
    stop("df_list must be a list of data frames")
  }
  
  # Extract unique genes from all data frames
  all_genes <- sort(unique(unlist(BiocParallel::bplapply(df_list, function(data) {
    unique(c(as.character(data[[1]]), as.character(data[[2]])))
  }, BPPARAM = BiocParallel::MulticoreParam(nCores)))))
  
  # Create a template adjacency matrix
  template_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_genes),
                            dimnames = list(all_genes, all_genes))
  
  # Process each data frame in parallel
  adjacency_matrix_list <- BiocParallel::bplapply(df_list, function(data) {
    adjacency_matrix <- template_matrix  # Copy template structure
    
    for (i in seq_len(nrow(data))) {
      gene1 <- as.character(data[i, 1])
      gene2 <- as.character(data[i, 2])
      weight <- as.numeric(data[i, 3])
      
      if (!is.na(weight) && gene1 %in% rownames(adjacency_matrix) && gene2 %in% colnames(adjacency_matrix)) {
        adjacency_matrix[gene1, gene2] <- weight
      }
    }
    
    diag(adjacency_matrix) <- 0  # Remove self-loops
    return(adjacency_matrix)
  }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  
  return(adjacency_matrix_list)
}

