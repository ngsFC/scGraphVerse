#' Edge Mining of Gene Interactions Using PubMed
#'
#' This function compares predicted gene adjacency matrices to a ground truth matrix,
#' identifies gene pairs that are predicted as interacting (value 1) but are absent (value 0)
#' in the ground truth, and then queries PubMed for literature evidence of interactions between
#' these gene pairs. For each pair, it returns the number of PubMed hits and the associated PubMed IDs.
#'
#' @param predicted_list A list of predicted gene adjacency matrices. Each matrix must have row names
#'   and column names representing gene names.
#' @param ground_truth A ground truth gene adjacency matrix (0/1) with row names and column names.
#' @param delay Numeric. The delay in seconds between consecutive PubMed queries to avoid overwhelming
#'   the NCBI servers. Default is 0.5.
#' @param query_field Character. The PubMed search field to be used in the query (e.g., "Title/Abstract").
#'   Default is "Title/Abstract".
#'
#' @return A list of data frames, one for each predicted matrix. Each data frame contains:
#'   \describe{
#'     \item{gene1}{The first gene in the gene pair.}
#'     \item{gene2}{The second gene in the gene pair.}
#'     \item{pubmed_hits}{The number of PubMed hits returned for the gene pair.}
#'     \item{PMIDs}{A comma-separated string of PubMed IDs associated with the hits (or \code{NA} if none).}
#'   }
#'
#' @details
#' The function identifies discrepant gene pairs where the predicted matrix indicates an interaction (1)
#' but the ground truth does not (0). For each discrepant pair, a query is constructed and sent to PubMed
#' using the \code{rentrez} package. The query first retrieves the count of hits and, if hits are found,
#' a second query is made to retrieve all PubMed IDs. A delay is introduced between queries to comply with
#' NCBI usage guidelines.
#'
#' @examples
#' \dontrun{
#' # Load predicted matrices and ground truth from CSV files (example):
#' predicted1 <- as.matrix(read.csv("predicted_matrix1.csv", row.names = 1))
#' predicted2 <- as.matrix(read.csv("predicted_matrix2.csv", row.names = 1))
#' ground_truth <- as.matrix(read.csv("ground_truth.csv", row.names = 1))
#'
#' # Combine predicted matrices into a list
#' predicted_list <- list(matrix1 = predicted1, matrix2 = predicted2)
#'
#' # Run the edge mining function
#' result <- edge_mining(predicted_list, ground_truth, delay = 0.5, query_field = "Title/Abstract")
#'
#' # View the result for the first matrix
#' print(result[[1]])
#' }
#'
#' @import rentrez
#' @export

edge_mining <- function(predicted_list, ground_truth, delay = 0.5, query_field = "Title/Abstract", 
                        query_edge_types = c("TP", "FP", "FN")) {
  
  # Check that the ground_truth matrix has gene names in its row and column names.
  if (is.null(rownames(ground_truth)) || is.null(colnames(ground_truth))) {
    stop("The ground_truth matrix must have row and column names representing gene names.")
  }
  
  # Initialize a list to store the results for each predicted matrix.
  results_list <- list()
  
  # Loop over each predicted matrix in the list.
  for (i in seq_along(predicted_list)) {
    
    predicted <- predicted_list[[i]]
    
    # Check that the predicted matrix has gene names.
    if (is.null(rownames(predicted)) || is.null(colnames(predicted))) {
      stop(paste("Predicted matrix at index", i, "does not have proper row and column names."))
    }
    
    # Identify all gene pairs that are present in either predicted or ground_truth.
    # (This excludes pairs where both matrices have a 0.)
    indices <- which((predicted == 1) | (ground_truth == 1), arr.ind = TRUE)
    
    # If no gene pairs are found, store an empty data frame.
    if (nrow(indices) == 0) {
      message("No gene pairs found for predicted matrix index ", i)
      results_list[[i]] <- data.frame(gene1 = character(0),
                                      gene2 = character(0),
                                      edge_type = character(0),
                                      pubmed_hits = integer(0),
                                      PMIDs = character(0),
                                      stringsAsFactors = FALSE)
      next
    }
    
    # Create a data frame of the gene pairs.
    gene_pairs <- data.frame(
      gene1 = rownames(predicted)[indices[, "row"]],
      gene2 = colnames(predicted)[indices[, "col"]],
      stringsAsFactors = FALSE
    )
    
    # Determine the edge type (TP, FP, or FN) for each gene pair.
    gene_pairs$edge_type <- NA_character_
    for (k in seq_len(nrow(gene_pairs))) {
      r <- indices[k, "row"]
      c <- indices[k, "col"]
      
      if (predicted[r, c] == 1 && ground_truth[r, c] == 1) {
        gene_pairs$edge_type[k] <- "TP"
      } else if (predicted[r, c] == 1 && ground_truth[r, c] == 0) {
        gene_pairs$edge_type[k] <- "FP"
      } else if (predicted[r, c] == 0 && ground_truth[r, c] == 1) {
        gene_pairs$edge_type[k] <- "FN"
      }
    }
    
    # Filter gene pairs based on the query_edge_types parameter.
    gene_pairs <- gene_pairs[gene_pairs$edge_type %in% query_edge_types, , drop = FALSE]
    
    # Initialize columns for PubMed query results.
    gene_pairs$pubmed_hits <- NA_integer_
    gene_pairs$PMIDs <- NA_character_
    
    # Loop over each gene pair and query PubMed.
    for (j in seq_len(nrow(gene_pairs))) {
      gene1 <- gene_pairs$gene1[j]
      gene2 <- gene_pairs$gene2[j]
      
      # Build the query string.
      query <- paste0(gene1, "[", query_field, "] AND ", gene2, "[", query_field, "]")
      message("Matrix ", i, ": Querying PubMed for: ", query)
      
      # Use tryCatch to safely query PubMed.
      result <- tryCatch({
        # First, do a query to get the count.
        search_res <- entrez_search(db = "pubmed", term = query, retmax = 0)
        hit_count <- as.numeric(search_res$count)
        
        # If hits are found, perform a second query to retrieve all PMIDs.
        if (hit_count > 0) {
          search_res <- entrez_search(db = "pubmed", term = query, retmax = hit_count)
        }
        search_res
      }, error = function(e) {
        message("Error querying PubMed for gene pair: ", gene1, " and ", gene2, ". Error: ", e$message)
        return(NULL)
      })
      
      # If the query was successful, store the hit count and PMIDs.
      if (!is.null(result)) {
        gene_pairs$pubmed_hits[j] <- as.numeric(result$count)
        if (length(result$ids) > 0) {
          # Collapse PMIDs into a comma-separated string.
          gene_pairs$PMIDs[j] <- paste(result$ids, collapse = ",")
        } else {
          gene_pairs$PMIDs[j] <- NA_character_
        }
      }
      
      # Pause briefly between queries to be polite to the NCBI servers.
      Sys.sleep(delay)
    }
    
    # Save the results for this matrix.
    results_list[[i]] <- gene_pairs
  }
  
  # If the input list has names, assign them to the results list.
  if (!is.null(names(predicted_list))) {
    names(results_list) <- names(predicted_list)
  }
  
  return(results_list)
}
