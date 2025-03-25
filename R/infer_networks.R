#' Infer Gene Regulatory Networks Using Various Methods
#'
#' This function infers gene regulatory networks (GRNs) from a list of expression matrices
#' using the selected method: "GENIE3", "GRNBoost2", "ZILGM", "JRF", or "PCzinb".
#'
#' @param count_matrices_list A list of expression matrices (genes × cells) or Seurat / SingleCellExperiment objects.
#' @param method Character. One of: "GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb".
#' @param adjm Optional adjacency matrix (used only for ZILGM and PCzinb row/col naming).
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
                           nCores = BiocParallel::bpworkers(BiocParallel::bpparam()) - 1,
                           grnboost_modules = NULL) {

  method <- match.arg(method, choices = c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"))

  # Normalize all inputs to matrices
  count_matrices_list <- lapply(count_matrices_list, function(obj) {
    if (inherits(obj, "Seurat")) {
      expr <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
      as.matrix(expr)
    } else if (inherits(obj, "SingleCellExperiment")) {
      expr <- SummarizedExperiment::assay(obj, "counts")
      as.matrix(expr)
    } else {
      as.matrix(obj)
    }
  })

  # --- GENIE3 ---
  if (method == "GENIE3") {
    return(BiocParallel::bplapply(count_matrices_list, function(mat) {
      adj <- GENIE3::GENIE3(mat, nCores = nCores)
      GENIE3::getLinkList(adj)
    }, BPPARAM = BiocParallel::MulticoreParam(nCores)))
  }

  # --- GRNBoost2 ---
  if (method == "GRNBoost2") {
    if (is.null(grnboost_modules)) {
      stop("For method 'GRNBoost2', please provide Python modules via `init_py()`.")
    }

    return(BiocParallel::bplapply(seq_along(count_matrices_list), function(i) {
      mat <- count_matrices_list[[i]]
      df <- as.data.frame(t(mat))
      genes <- colnames(df)
      df_pandas <- grnboost_modules$pandas$DataFrame(
        data = as.matrix(df),
        columns = genes,
        index = rownames(df)
      )
      grnboost_modules$arboreto$grnboost2(df_pandas, gene_names = genes)
    }, BPPARAM = BiocParallel::MulticoreParam(nCores)))
  }

  # --- ZILGM ---
  if (method == "ZILGM") {
    zilgm_fits <- BiocParallel::bplapply(count_matrices_list, function(mat) {
      lambda_max <- ZILGM::find_lammax(t(mat))
      lambda_seq <- exp(seq(log(lambda_max), log(1e-4 * lambda_max), length.out = 50))
      fit <- ZILGM::zilgm(
        X = t(mat),
        lambda = lambda_seq,
        family = "NBII",
        update_type = "IRLS",
        do_boot = TRUE,
        boot_num = 10,
        sym = "OR",
        nCores = nCores
      )
      fit
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))

    network_results <- lapply(zilgm_fits, function(fit) {
      mat <- fit$network[[fit$opt_index]]
      if (!is.null(adjm)) {
        rownames(mat) <- rownames(adjm)
        colnames(mat) <- colnames(adjm)
      }
      mat
    })

    lambda_results <- lapply(zilgm_fits, function(fit) {
      nets <- lapply(fit$network, as.matrix)
      names(nets) <- fit$lambda
      if (!is.null(adjm)) {
        for (i in seq_along(nets)) {
          rownames(nets[[i]]) <- rownames(adjm)
          colnames(nets[[i]]) <- colnames(adjm)
        }
      }
      nets
    })

    return(list(network_results = network_results, lambda_results = lambda_results))
  }

  # --- PCzinb ---
  if (method == "PCzinb") {
    return(BiocParallel::bplapply(count_matrices_list, function(mat) {
      net <- learn2count::PCzinb(t(mat), method = "zinb1", maxcard = 2)
      if (!is.null(adjm)) {
        rownames(net) <- rownames(adjm)
        colnames(net) <- colnames(adjm)
      }
      net
    }, BPPARAM = BiocParallel::MulticoreParam(nCores)))
  }

  # --- JRF (joint inference across all matrices) ---
  if (method == "JRF") {
    norm_list <- BiocParallel::bplapply(count_matrices_list, function(x) {
      (x - mean(x)) / sd(x)
    }, BPPARAM = BiocParallel::MulticoreParam(nCores))

    rf <- JRF::JRF(
      X = norm_list,
      genes.name = rownames(norm_list[[1]]),
      ntree = 500,
      mtry = round(sqrt(nrow(norm_list[[1]]) - 1))
    )
    return(list(rf))
  }
}

