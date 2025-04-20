#' Infer Gene Regulatory Networks Using Various Methods
#'
#' This function infers gene regulatory networks (GRNs) from a list of expression matrices
#' using the selected method: "GENIE3", "GRNBoost2", "ZILGM", "JRF", or "PCzinb".
#'
#' @param count_matrices_list A list of expression matrices (genes × cells) or Seurat / SingleCellExperiment objects.
#' @param method Character. One of: "GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb".
#' @param adjm Optional adjacency matrix (used to rename output matrices).
#' @param nCores Number of CPU cores to use for parallel computation.
#' @param grnboost_modules Python modules list returned by `init_py()`, required for GRNBoost2.
#'
#' @return A list of inferred networks (format varies by method).
#'
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
#' @export
infer_networks <- function(count_matrices_list,
                           method = "GENIE3",
                           adjm = NULL,
                           total_cores = BiocParallel::bpworkers(BiocParallel::bpparam()),
                           grnboost_modules = NULL,
                           seed = 123) {
  
  method <- match.arg(method, c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"))
  
  # Define outer and inner cores based on method
  n_matrices <- length(count_matrices_list)
  if (method == "GRNBoost2") {
    param_outer <- BiocParallel::SerialParam()
    nCores_inner <- total_cores
  } else {
    nCores_outer <- min(total_cores, n_matrices)
    nCores_inner <- max(floor(total_cores / nCores_outer), 1)
    param_outer <- BiocParallel::MulticoreParam(workers = nCores_outer, RNGseed = seed)
  }
  
  count_matrices_list <- lapply(count_matrices_list, function(obj) {
    if (inherits(obj, "Seurat")) {
      as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
    } else if (inherits(obj, "SingleCellExperiment")) {
      as.matrix(SummarizedExperiment::assay(obj, "counts"))
    } else {
      as.matrix(obj)
    }
  })
  
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
      
      df_pandas <- grnboost_modules$pandas$DataFrame(
        data = as.matrix(df),
        columns = genes,
        index = rownames(df)
      )
      
      result_py <- grnboost_modules$arboreto$grnboost2(
        expression_data = df_pandas,
        gene_names = genes
      )
      
      result_r <- reticulate::py_to_r(result_py)
      if (is.data.frame(result_r)) {
        rownames(result_r) <- NULL
      }
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
      
    } else if (method == "JRF") {
      param_inner <- BiocParallel::MulticoreParam(workers = nCores_inner, RNGseed = task_seed)
      #norm_mat <- (mat - mean(mat)) / sd(mat)
      #mat <- log(mat+1)
      norm_mat <- t(scale(t(mat)))
      #library(edgeR)
      #dge   <- DGEList(counts=mat)
      #dge   <- calcNormFactors(dge, method="TMM")
      #norm_mat <- cpm(dge, log=TRUE, prior.count=1)
      rf <- JRF::JRF(X = list(norm_mat),
                     genes.name = rownames(norm_mat),
                     ntree = 500,
                     mtry = round(sqrt(nrow(norm_mat) - 1)))
      rf
    }
  }, BPPARAM = param_outer)
}
