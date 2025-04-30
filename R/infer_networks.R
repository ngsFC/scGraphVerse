#' Infer Gene Regulatory Networks from Expression Matrices
#'
#' Infers weighted gene regulatory networks (GRNs) from one or more expression matrices
#' using different inference methods: \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"}, \code{"JRF"}, or \code{"PCzinb"}.
#'
#' @param count_matrices_list A list of expression matrices (genes × cells) or
#'   \linkS4class{Seurat} or \linkS4class{SingleCellExperiment} objects.
#' @param method Character string. Inference method to use. One of:
#'   \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"}, \code{"JRF"}, or \code{"PCzinb"}.
#' @param adjm Optional. Reference adjacency matrix for matching dimensions when using \code{"ZILGM"} or \code{"PCzinb"}.
#' @param nCores Integer. Number of CPU cores to use for parallelization.
#'   Defaults to the number of workers in the current \pkg{BiocParallel} backend.
#' @param grnboost_modules Python modules required for \code{GRNBoost2} (created via \pkg{reticulate}).
#'
#' @return
#' A list of inferred networks:
#' \itemize{
#'   \item For \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"}, and \code{"PCzinb"}, a list of inferred network objects (edge lists or adjacency matrices).
#'   \item For \code{"JRF"}, a list of data frames with inferred edge lists for each condition or dataset.
#' }
#'
#' @details
#' Each expression matrix is preprocessed automatically depending on its object type (\code{Seurat},
#' \code{SingleCellExperiment}, or plain matrix).
#'
#' Parallelization behavior:
#' \itemize{
#'   \item \strong{GENIE3} and \strong{ZILGM}: No external parallelization; internal \code{nCores} parameter controls computation.
#'   \item \strong{GRNBoost2} and \strong{PCzinb}: Parallelized across matrices using \pkg{BiocParallel}.
#'   \item \strong{JRF}: Joint modeling of all matrices together; internal parallelization across random forest trees using \pkg{doParallel}.
#' }
#'
#' Methods are based on:
#' \itemize{
#'   \item \strong{GENIE3}: Random Forest-based inference (Huynh-Thu et al., 2010).
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
#' set.seed(123)
#'
#' # Simulate two small expression matrices
#' mat1 <- matrix(rpois(100, lambda = 5), nrow = 10)
#' mat2 <- matrix(rpois(100, lambda = 5), nrow = 10)
#' rownames(mat1) <- paste0("Gene", 1:10)
#' rownames(mat2) <- paste0("Gene", 1:10)
#'
#' # Infer networks using GENIE3
#' networks <- infer_networks(
#'   count_matrices_list = list(mat1, mat2),
#'   method = "GENIE3",
#'   nCores = 2
#' )
#'
#' # Inspect first inferred network
#' head(networks[[1]])
infer_networks <- function(count_matrices_list,
                           method = "GENIE3",
                           adjm = NULL,
                           nCores = BiocParallel::bpworkers(BiocParallel::bpparam()),
                           grnboost_modules = NULL) {
  method <- match.arg(method, c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"))

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

  # GENIE3 and ZILGM: serial loop
  if (method %in% c("GENIE3", "ZILGM")) {
    results <- vector("list", length = n_matrices)
    for (i in seq_len(n_matrices)) {
      mat <- count_matrices_list[[i]]

      if (method == "GENIE3") {
        adj <- GENIE3::GENIE3(mat, nCores = nCores)
        results[[i]] <- GENIE3::getLinkList(adj)
      } else if (method == "ZILGM") {
        lambda_max <- ZILGM::find_lammax(t(mat))
        lambda_seq <- exp(seq(log(lambda_max), log(1e-4 * lambda_max), length.out = 50))
        fit <- ZILGM::zilgm(
          X = t(mat), lambda = lambda_seq, family = "NBII",
          update_type = "IRLS", do_boot = TRUE, boot_num = 10,
          sym = "OR", nCores = nCores
        )
        adj <- fit$network[[fit$opt_index]]
        dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
        results[[i]] <- adj
      }
    }
    return(results)
  }

  # JRF: all matrices together, parallel inside
  if (method == "JRF") {
    norm_list <- lapply(count_matrices_list, function(mat) {
      t(scale(t(mat)))
    })

    clust <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(clust), add = TRUE)
    doParallel::registerDoParallel(clust)

    rf <- JRF::JRF(
      X = norm_list,
      genes.name = rownames(norm_list[[1]]),
      ntree = 500,
      mtry = round(sqrt(nrow(norm_list[[1]]) - 1))
    )

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

  # GRNBoost2 and PCzinb: parallel across matrices
  param_outer <- BiocParallel::MulticoreParam(workers = nCores)
  BiocParallel::bplapply(seq_along(count_matrices_list), function(i) {
    mat <- count_matrices_list[[i]]

    if (method == "GRNBoost2") {
      if (is.null(grnboost_modules)) stop("Provide grnboost_modules for GRNBoost2.")

      df <- as.data.frame(t(mat))
      genes <- colnames(df)
      rownames(df) <- make.unique(rownames(df))
      df_pandas <- grnboost_modules$pandas$DataFrame(data = as.matrix(df), columns = genes, index = rownames(df))

      result_py <- grnboost_modules$arboreto$grnboost2(expression_data = df_pandas, gene_names = genes)
      result_r <- reticulate::py_to_r(result_py)
      if (is.data.frame(result_r)) rownames(result_r) <- NULL
      result_r
    } else if (method == "PCzinb") {
      adj <- learn2count::PCzinb(t(mat), method = "zinb1", maxcard = 2)
      dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
      adj
    }
  }, BPPARAM = param_outer)
}
