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
#' @param nCores The number of cores to use for parallel computation. Default is the number of available cores.
#'
#' @return A list of data frames, one for each predicted matrix. Each data frame contains:
#'   \describe{
#'     \item{gene1}{The first gene in the gene pair.}
#'     \item{gene2}{The second gene in the gene pair.}
#'     \item{pubmed_hits}{The number of PubMed hits returned for the gene pair.}
#'     \item{PMIDs}{A comma-separated string of PubMed IDs associated with the hits (or \code{NA} if none).}
#'   }
#'
#' @import BiocParallel rentrez
#' @export

library(BiocParallel)
library(rentrez)

edge_mining <- function(predicted_list, ground_truth, delay = 0.5, query_field = "Title/Abstract", 
                        query_edge_types = c("TP", "FP", "FN"),
                        nCores = BiocParallel::bpworkers(BiocParallel::bpparam())) {
  
  if (!is.list(predicted_list) || !is.matrix(ground_truth)) {
    stop("predicted_list must be a list of matrices and ground_truth must be a matrix")
  }
  
  results_list <- BiocParallel::bplapply(seq_along(predicted_list), function(i) {
    predicted <- predicted_list[[i]]
    
    if (!is.matrix(predicted) || is.null(rownames(predicted)) || is.null(colnames(predicted))) {
      stop(paste("Predicted matrix at index", i, "does not have proper row and column names."))
    }
    
    indices <- which(((predicted == 1) | (ground_truth == 1)) & upper.tri(predicted), arr.ind = TRUE)
    if (nrow(indices) == 0) {
      return(data.frame(gene1 = character(0), gene2 = character(0), edge_type = character(0),
                        pubmed_hits = integer(0), PMIDs = character(0), stringsAsFactors = FALSE))
    }
    
    gene_pairs <- data.frame(
      gene1 = rownames(predicted)[indices[, "row"]],
      gene2 = colnames(predicted)[indices[, "col"]],
      stringsAsFactors = FALSE
    )
    
    gene_pairs$edge_type <- ifelse(predicted[indices] == 1 & ground_truth[indices] == 1, "TP",
                             ifelse(predicted[indices] == 1 & ground_truth[indices] == 0, "FP",
                             "FN"))
    
    gene_pairs <- gene_pairs[gene_pairs$edge_type %in% query_edge_types, , drop = FALSE]
    
    if (nrow(gene_pairs) == 0) return(gene_pairs)
    
    pubmed_results <- BiocParallel::bplapply(seq_len(nrow(gene_pairs)), function(j) {
      gene1 <- gene_pairs$gene1[j]
      gene2 <- gene_pairs$gene2[j]
      query <- paste0(gene1, "[", query_field, "] AND ", gene2, "[", query_field, "]")
      
      result <- tryCatch({
        search_res <- entrez_search(db = "pubmed", term = query, retmax = 0)
        hit_count <- as.numeric(search_res$count)
        
        if (hit_count > 0) {
          search_res <- entrez_search(db = "pubmed", term = query, retmax = hit_count)
        }
        search_res
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(result)) {
        pubmed_hits <- as.numeric(result$count)
        pmids <- if (length(result$ids) > 0) paste(result$ids, collapse = ",") else NA_character_
      } else {
        pubmed_hits <- NA_integer_
        pmids <- NA_character_
      }
      
      Sys.sleep(delay)  # Respect NCBI guidelines
      return(data.frame(pubmed_hits, PMIDs = pmids))
    }, BPPARAM = MulticoreParam(nCores))
    
    pubmed_results <- do.call(rbind, pubmed_results)
    gene_pairs$pubmed_hits <- pubmed_results$pubmed_hits
    gene_pairs$PMIDs <- pubmed_results$PMIDs
    
    return(gene_pairs)
  }, BPPARAM = MulticoreParam(nCores))
  
  if (!is.null(names(predicted_list))) {
    names(results_list) <- names(predicted_list)
  }
  
  return(results_list)
}

