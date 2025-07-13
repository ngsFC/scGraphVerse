#' Infer Gene Regulatory Networks from Expression Matrices
#'
#' Infers weighted gene regulatory networks (GRNs) from one or more
#' expression matrices using different inference methods:
#' \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"},
#' \code{"JRF"}, or \code{"PCzinb"}.
#'
#' @param count_matrices_list A list of expression matrices (genes × cells)
#'   or \linkS4class{Seurat} or \linkS4class{SingleCellExperiment}
#'   objects.
#' @param method Character string. Inference method to use. One of:
#'   \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"},
#'   \code{"JRF"}, or \code{"PCzinb"}.
#' @param adjm Optional. Reference adjacency matrix for matching dimensions
#'   when using \code{"ZILGM"} or \code{"PCzinb"}.
#' @param nCores Integer. Number of CPU cores to use for
#'   parallelization. Defaults to the number of workers in the current
#'   \pkg{BiocParallel} backend.
#' @param grnboost_modules Python modules required for \code{GRNBoost2}
#'   (created via \pkg{reticulate}).
#' @param genie3_params List of parameters for GENIE3 method:
#'   \itemize{
#'     \item \code{regulators}: Vector of regulator gene names (default: all)
#'     \item \code{targets}: Vector of target gene names (default: all genes)
#'     \item \code{treeMethod}: "RF" or "ET" (default: "RF")
#'     \item \code{K}: Number of candidate regulators (default: "sqrt")
#'     \item \code{nTrees}: Number of trees per ensemble (default: 1000)
#'     \item \code{seed}: Random seed for reproducibility (default: NULL)
#'   }
#' @param grnboost2_params List of parameters for GRNBoost2 method:
#'   \itemize{
#'     \item \code{tf_names}: Vector of transcription factor names (default:all)
#'     \item \code{gene_names}: Vector of target gene names (default: all)
#'     \item \code{client_or_address}: Dask client or address (default: NULL)
#'     \item \code{seed}: Random seed for reproducibility (default: NULL)
#'   }
#' @param zilgm_params List of parameters for ZILGM method:
#'   \itemize{
#'     \item \code{lambda}: Regularization parameter (default: 0.1)
#'     \item \code{alpha}: Elastic net mixing parameter (default: 1)
#'     \item \code{max_iter}: Maximum iterations (default: 100)
#'     \item \code{tol}: Convergence tolerance (default: 1e-4)
#'   }
#' @param jrf_params List of parameters for JRF method:
#'   \itemize{
#'     \item \code{ntree}: Number of trees (default: 500)
#'     \item \code{mtry}: Number of variables to sample at each split 
#'     (default: sqrt(p))
#'     \item \code{nodesize}: Minimum node size (default: 5)
#'     \item \code{maxnodes}: Maximum number of nodes (default: NULL)
#'   }
#' @param pczinb_params List of parameters for PCzinb method:
#'   \itemize{
#'     \item \code{gamma}: Regularization parameter (default: 0.1)
#'     \item \code{beta}: Beta parameter (default: 0.1)
#'     \item \code{max_iter}: Maximum iterations (default: 100)
#'     \item \code{tol}: Convergence tolerance (default: 1e-4)
#'   }
#' @param verbose Logical. If TRUE, display progress messages. Default: FALSE.
#' @param seed Integer. Random seed for reproducibility. Default: NULL.
#'
#' @return A list of inferred networks:
#'   \itemize{
#'     \item For \code{"GENIE3"}, \code{"GRNBoost2"}, \code{"ZILGM"},
#'       and \code{"PCzinb"}, a list of inferred network objects (edge
#'       lists or adjacency matrices).
#'     \item For \code{"JRF"}, a list of data frames with inferred edge
#'       lists for each condition or dataset.
#'   }
#'
#' @details Each expression matrix is preprocessed automatically depending
#'   on its object type (\code{Seurat}, \code{SingleCellExperiment}, or
#'   plain matrix).
#'
#'   Parallelization behavior:
#'   \itemize{
#'     \item \strong{GENIE3} and \strong{ZILGM}: No external
#'       parallelization; internal \code{nCores} parameter controls
#'       computation.
#'     \item \strong{GRNBoost2} and \strong{PCzinb}: Parallelized across
#'       matrices using \pkg{BiocParallel}.
#'     \item \strong{JRF}: Joint modeling of all matrices together;
#'       internal parallelization across random forest trees using
#'       \pkg{doParallel}.
#'   }
#'
#'   Methods are based on:
#'   \itemize{
#'     \item \strong{GENIE3}: Random Forest-based inference (Huynh-Thu et
#'       al., 2010).
#'     \item \strong{GRNBoost2}: Gradient boosting trees using arboreto
#'       (Moerman et al., 2019).
#'     \item \strong{ZILGM}: Zero-Inflated Graphical Models for scRNA-seq
#'       (Zhang et al., 2021).
#'     \item \strong{JRF}: Joint Random Forests across multiple conditions
#'       (Petralia et al., 2015).
#'     \item \strong{PCzinb}: Pairwise correlation under ZINB models
#'       (Nguyen et al., 2023).
#'   }
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam bpworkers
#'   bpparam
#' @importFrom Seurat GetAssayData
#' @importFrom SummarizedExperiment assay
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @examples
#' data("count_matrices")
#'
#' networks <- infer_networks(
#'     count_matrices_list = count_matrices,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
infer_networks <- function(
    count_matrices_list,
    method = c("GENIE3", "GRNBoost2", "ZILGM", "JRF", "PCzinb"),
    adjm = NULL,
    nCores = 1,
    grnboost_modules = NULL,
    genie3_params = list(),
    grnboost2_params = list(),
    zilgm_params = list(),
    jrf_params = list(),
    pczinb_params = list(),
    verbose = FALSE,
    seed = NULL) {
    method <- match.arg(method)
    count_matrices_list <- .convert_counts_list(count_matrices_list)
    n_matrices <- length(count_matrices_list)
    
    # Note: Seed handling is method-specific and implemented within each method
    # to comply with Bioconductor guidelines
    
    # Merge method-specific parameters with defaults
    genie3_params <- .merge_genie3_params(genie3_params)
    grnboost2_params <- .merge_grnboost2_params(grnboost2_params)
    zilgm_params <- .merge_zilgm_params(zilgm_params)
    jrf_params <- .merge_jrf_params(jrf_params)
    pczinb_params <- .merge_pczinb_params(pczinb_params)

    if (method %in% c("GENIE3", "ZILGM")) {
        results <- vector("list", n_matrices)
        for (i in seq_len(n_matrices)) {
            mat <- count_matrices_list[[i]]
            if (method == "GENIE3") {
                if (verbose) {
                    message("Running GENIE3 on matrix ", i, "/", n_matrices)
                }
                results[[i]] <- .run_genie3(mat, nCores, genie3_params)
            } else {
                if (verbose) {
                    message("Running ZILGM on matrix ", i, "/", n_matrices)
                }
                results[[i]] <- .run_zilgm(mat, adjm, nCores, zilgm_params)
            }
        }
        return(results)
    }

    if (method == "JRF") {
        norm_list <- lapply(
            count_matrices_list,
            function(mat) t(scale(t(mat)))
        )
        if (verbose) message("Running JRF on all matrices jointly")
        return(.run_jrf(norm_list, nCores, jrf_params))
    }

    if (method == "GRNBoost2") {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop(
                "'reticulate' package is required for method = 'GRNBoost2'.\n",
                "Install with: install.packages('reticulate')",
                call. = FALSE
            )
        }
        
        # Check if Python and arboreto are available
        python_available <- tryCatch({
            reticulate::py_available(initialize = TRUE)
        }, error = function(e) FALSE)
        
        if (!python_available) {
            stop(
                "Python is not available or not properly configured.\n",
                "Please ensure Python is installed and accessible.\n",
                "You may need to restart R after installing Python.\n",
                "For setup help, see: vignette('python-setup', package = 'scGraphVerse')",
                call. = FALSE
            )
        }
        
        arboreto_available <- tryCatch({
            reticulate::py_module_available("arboreto")
        }, error = function(e) FALSE)
        
        if (!arboreto_available) {
            stop(
                "Python package 'arboreto' is required for GRNBoost2.\n",
                "Install options:\n",
                "  1. Automatic: setup_grnboost2()\n",
                "  2. Manual pip: pip install arboreto\n",
                "  3. Manual conda: conda install -c bioconda arboreto\n",
                "  4. From R: reticulate::py_install('arboreto')\n",
                "For detailed setup instructions, see: vignette('python-setup', package = 'scGraphVerse')",
                call. = FALSE
            )
        }
        if (verbose) message("Running GRNBoost2 on ", n_matrices, " matrices")
        return(.run_parallel_networks(
            count_matrices_list,
            method,
            nCores,
            adjm,
            grnboost_modules,
            grnboost2_params
        ))
    }

    if (method == "PCzinb") {
        if (!requireNamespace("learn2count", quietly = TRUE)) {
            stop(
                "Package 'learn2count' is required for method = 'PCzinb'.\n",
                "Please install it: BiocManager::install('drisso/learn2count')",
                call. = FALSE
            )
        }
        if (verbose) message("Running PCzinb on ", n_matrices, " matrices")
        return(.run_parallel_networks(
            count_matrices_list,
            method,
            nCores,
            adjm,
            grnboost_modules,
            pczinb_params
        ))
    }
}

