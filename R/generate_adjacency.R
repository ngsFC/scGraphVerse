#' Generate Adjacency Matrices from a List of Data Frames
#'
#' This function generates adjacency matrices based on the information provided in 
#' a list of data frames. Each data frame is expected to have three columns: 
#' the first and second columns represent gene pairs, and the third column represents 
#' the weight of the interaction between the genes. The adjacency matrices are constructed 
#' based on these gene pairs and weights, matching the row and column names of the 
#' provided `ground.truth` adjacency matrix.
#'
#' @param df_list A list of data frames, each containing three columns: 
#'   - First column: Gene1 (character)
#'   - Second column: Gene2 (character)
#'   - Third column: Weight (numeric)
#' @param ground.truth An adjacency matrix (or data frame) with the ground truth 
#'   values, used for the dimensions and gene names of the resulting adjacency matrices.
#' @return A list of adjacency matrices, one for each data frame in `df_list`.
#' @details 
#' The function creates an adjacency matrix for each data frame in `df_list` by placing 
#' the interaction weights in the correct positions based on the gene pairings. The 
#' diagonal of the resulting matrices is set to zero to ensure no self-interactions are included.
#' 
#' @examples
#' \dontrun{
#' # Example data frame
#' df1 <- data.frame(Gene1 = c("GeneA", "GeneB"), Gene2 = c("GeneB", "GeneC"), Weight = c(0.5, 0.8))
#' df_list <- list(df1)
#' ground_truth <- matrix(0, nrow = 3, ncol = 3, dimnames = list(c("GeneA", "GeneB", "GeneC"), c("GeneA", "GeneB", "GeneC")))
#' adjacency_matrices <- generate_adjacency(df_list, ground_truth)
#' }
#' @export
generate_adjacency <- function(df_list, ground.truth = NULL) {
  adjacency_matrix_list <- list()
  
  # If ground.truth is not supplied, derive gene names from df_list
  if (is.null(ground.truth)) {
    # Extract gene names from the first two columns of all data frames
    all_genes <- unique(unlist(lapply(df_list, function(data) {
      unique(c(as.character(data[, 1]), as.character(data[, 2])))
    })))
    # Sort genes in alphabetical order
    all_genes <- sort(all_genes)
    
    # Create a template matrix with rows and columns in alphabetical order
    template_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_genes))
    rownames(template_matrix) <- all_genes
    colnames(template_matrix) <- all_genes
  } else {
    # Use the supplied ground.truth as the template
    template_matrix <- ground.truth
  }
  
  # Process each data frame in the list
  for (k in seq_along(df_list)) {
    data <- df_list[[k]]
    # Create an empty adjacency matrix with the same dimensions and names as template_matrix
    adjacency_matrix <- matrix(0, nrow = nrow(template_matrix), ncol = ncol(template_matrix))
    rownames(adjacency_matrix) <- rownames(template_matrix)
    colnames(adjacency_matrix) <- colnames(template_matrix)
    
    # Fill the matrix with weights from the current data frame
    for (i in 1:nrow(data)) {
      gene1 <- as.character(data[i, 1])
      gene2 <- as.character(data[i, 2])
      weight <- data[i, 3]
      if (gene1 %in% rownames(adjacency_matrix) && gene2 %in% colnames(adjacency_matrix)) {
        adjacency_matrix[gene1, gene2] <- weight
      }
    }
    # Remove self-loops by setting diagonal to 0
    diag(adjacency_matrix) <- 0
    
    # Add the created adjacency matrix to the list
    adjacency_matrix_list[[k]] <- adjacency_matrix
  }
  
  return(adjacency_matrix_list)
}
