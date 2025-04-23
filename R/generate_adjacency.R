#' Generate Adjacency Matrices from Gene Interaction Tables
#'
#' This function generates adjacency matrices from a list of data frames, where each 
#' data frame represents weighted gene interactions. Each row specifies an interaction 
#' between two genes, along with an associated weight.
#'
#' @param df_list A list of data frames, each containing three columns:
#'   \describe{
#'     \item{Gene1}{Character. First gene in the interaction.}
#'     \item{Gene2}{Character. Second gene in the interaction.}
#'     \item{Weight}{Numeric. Strength of the interaction from Gene1 to Gene2.}
#'   }
#' @param nCores Integer (optional). Number of cores to use for parallel processing. 
#'   Defaults to the number of available workers on the current BiocParallel backend.
#'
#' @return A list of square numeric matrices (adjacency matrices). Each matrix has 
#'   genes as row and column names, and interaction weights as values. Diagonal 
#'   entries are set to zero.
#'
#' @details 
#' This function aggregates all unique gene names across all input data frames to 
#' define the full matrix dimensions. For each input data frame, a corresponding 
#' adjacency matrix is constructed by placing weights at the appropriate gene pair 
#' coordinates. Diagonal values are zeroed to exclude self-interactions. Parallel 
#' processing is used via the \pkg{BiocParallel} framework for efficient computation.
#'
#' @examples
#' df1 <- data.frame(Gene1 = c("GeneA", "GeneB"),
#'                   Gene2 = c("GeneB", "GeneC"),
#'                   Weight = c(0.5, 0.8))
#' df2 <- data.frame(Gene1 = c("GeneC"), Gene2 = c("GeneA"), Weight = 1.0)
#' adjacency_list <- generate_adjacency(list(df1, df2))
#' 
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
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
    adjacency_matrix <- template_matrix 
    
    for (i in seq_len(nrow(data))) {
      gene1 <- as.character(data[i, 1])
      gene2 <- as.character(data[i, 2])
      weight <- as.numeric(data[i, 3])
      
      if (!is.na(weight) && gene1 %in% rownames(adjacency_matrix) && gene2 %in% colnames(adjacency_matrix)) {
        adjacency_matrix[gene1, gene2] <- weight
      }
    }
    
    diag(adjacency_matrix) <- 0 
    return(adjacency_matrix)
  }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  
  return(adjacency_matrix_list)
}

