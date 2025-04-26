#' Infer Gene Regulatory Networks from Expression Matrices
#'
#' Infers weighted gene regulatory networks from one or more expression matrices
#' using different inference methods: \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"}, \code{"JRF"}, or \code{"PCzinb"}.
#'
#' @param count_matrices_list A list of expression matrices (genes × cells) or 
#'   \linkS4class{Seurat} or \linkS4class{SingleCellExperiment} objects.
#' @param method Character string. Inference method to use. One of:
#'   \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"}, \code{"JRF"}, or \code{"PCzinb"}.
#' @param adjm Optional. Reference adjacency matrix for matching dimensions when using \code{"ZILGM"} or \code{"PCzinb"}.
#' @param total_cores Integer. Total number of CPU cores for parallelization. 
#'   Defaults to the number of workers in the current \pkg{BiocParallel} backend.
#' @param grnboost_modules Python modules required for \code{GRNBoost2} (created via \pkg{reticulate}).
#' @param seed Integer. Random seed for reproducibility. Default is \code{123}.
#'
#' @return
#' \itemize{
#'   \item For \code{"JRF"}, a list of data frames with inferred edge lists for each dataset.
#'   \item For other methods, a list of inferred network objects (link lists or adjacency matrices).
#' }
#'
#' @details
#' Each expression matrix is preprocessed automatically depending on its object type (\code{Seurat}, 
#' \code{SingleCellExperiment}, or plain matrix).
#'
#' Parallelization across matrices is handled with \pkg{BiocParallel}. For \code{GENIE3} and \code{JRF},
#' nested parallelism is used: outer jobs distribute matrices and inner jobs distribute tree-based computations.
#'
#' The methods are based on:
#' \itemize{
#'   \item \strong{GENIE3}: Random Forest based inference (Huynh-Thu et al., 2010).
#'   \item \strong{GRNBoost2}: Gradient boosting trees using arboreto (Moerman et al., 2019).
#'   \item \strong{ZILGM}: Zero-Inflated Graphical Models for scRNA-seq (Zhang et al., 2021).
#'   \item \strong{JRF}: Joint Random Forests across multiple conditions (Petralia et al., 2015).
#'   \item \strong{PCzinb}: Pairwise correlation under ZINB models (Li et al., 2020).
#' }
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam bpworkers bpparam
#' @importFrom Seurat GetAssayData
#' @importFrom SummarizedExperiment assay
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @examples
#' set.seed(42)
#' # Simulate two small expression matrices
#' mat1 <- matrix(rpois(100, 5), nrow = 10)
#' mat2 <- matrix(rpois(100, 5), nrow = 10)
#' rownames(mat1) <- paste0("Gene", 1:10)
#' rownames(mat2) <- paste0("Gene", 1:10)
#'
#' # Infer networks using GENIE3
#' networks <- infer_networks(
#'   count_matrices_list = list(mat1, mat2),
#'   method = "GENIE3",
#'   total_cores = 2,
#'   seed = 123
#' )
#'
#' # Inspect first network
#' networks[[1]]

infer_networks <- function(count_matrices_list,
                           method = "GENIE3",
                           adjm = NULL,
                           total_cores = BiocParallel::bpworkers(BiocParallel::bpparam()),
                           grnboost_modules = NULL,
                           seed = 1234) {
  
  method <- match.arg(method, c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"))
  
  # Preprocess matrices
  count_matrices_list <- lapply(count_matrices_list, function(obj) {
    if (inherits(obj, "Seurat")) {
      as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
    } else if (inherits(obj, "SingleCellExperiment")) {
      as.matrix(SummarizedExperiment::assay(obj, "counts"))
    } else {
      as.matrix(obj)
    }
  })
  
  n_matrices <- length(count_matrices_list)
  nCores_outer <- min(total_cores, n_matrices)
  nCores_inner <- max(floor(total_cores / nCores_outer), 1)
  param_outer <- if (method == "GRNBoost2") BiocParallel::SerialParam() else BiocParallel::MulticoreParam(workers = nCores_outer, RNGseed = seed)
  
  if (method == "JRF") {
    norm_list <- lapply(count_matrices_list, function(mat) {
      t(scale(t(mat)))
    })
    
    param_inner <- BiocParallel::MulticoreParam(workers = nCores_inner, RNGseed = seed)
    clust <- parallel::makeCluster(nCores_inner)
    on.exit(parallel::stopCluster(clust), add = TRUE)
    doParallel::registerDoParallel(clust)
    
    # Run JRF joint network inference
    rf <- JRF::JRF(
      X = norm_list,
      genes.name = rownames(norm_list[[1]]),
      ntree = 500,
      mtry = round(sqrt(nrow(norm_list[[1]]) - 1))
    )
    
    # Split the result into a list of edge lists
    jrf_mat <- list(rf)
    importance_columns <- grep("importance", names(jrf_mat[[1]]), value = TRUE)
    
    jrf_list <- vector("list", length(importance_columns))
    for (i in seq_along(importance_columns)) {
      df <- jrf_mat[[1]][, c("gene1", "gene2", importance_columns[i])]
      names(df)[3] <- importance_columns[i]
      jrf_list[[i]] <- df
    }
    
    return(jrf_list)
  }
  
  BiocParallel::bplapply(seq_along(count_matrices_list), function(i) {
    mat <- count_matrices_list[[i]]
    task_seed <- as.integer(round(seed * 100 + i))
    set.seed(task_seed)
    
    if (method == "GENIE3") {
      adj <- GENIE3::GENIE3(mat, nCores = nCores_inner)
      GENIE3::getLinkList(adj)
      
    } else if (method == "GRNBoost2") {
      if (is.null(grnboost_modules)) stop("Provide grnboost_modules for GRNBoost2.")
      grnboost_modules$numpy$random$seed(as.integer(task_seed))
      
      df <- as.data.frame(t(mat))
      genes <- colnames(df)
      rownames(df) <- make.unique(rownames(df))
      df_pandas <- grnboost_modules$pandas$DataFrame(data = as.matrix(df), columns = genes, index = rownames(df))
      
      result_py <- grnboost_modules$arboreto$grnboost2(expression_data = df_pandas, gene_names = genes)
      result_r <- reticulate::py_to_r(result_py)
      if (is.data.frame(result_r)) rownames(result_r) <- NULL
      result_r
      
    } else if (method == "ZILGM") {
      lambda_max <- ZILGM::find_lammax(t(mat))
      lambda_seq <- exp(seq(log(lambda_max), log(1e-4 * lambda_max), length.out = 50))
      fit <- ZILGM::zilgm(X = t(mat), lambda = lambda_seq, family = "NBII",
                          update_type = "IRLS", do_boot = TRUE, boot_num = 10,
                          sym = "OR", nCores = nCores_inner)
      adj <- fit$network[[fit$opt_index]]
      dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
      adj
      
    } else if (method == "PCzinb") {
      param_inner <- BiocParallel::MulticoreParam(workers = nCores_inner, RNGseed = task_seed)
      adj <- learn2count::PCzinb(t(mat), method = "zinb1", maxcard = 2)
      dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
      adj
    }
  }, BPPARAM = param_outer)
}
