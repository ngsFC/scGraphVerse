infer_networks <- function(count_matrices_list,
                           method = "GENIE3",
                           adjm = NULL,
                           total_cores = BiocParallel::bpworkers(BiocParallel::bpparam()),
                           grnboost_modules = NULL,
                           seed = 123) {
  
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
  
  # JRF: joint modeling outside the loop
  if (method == "JRF") {
    norm_list <- lapply(count_matrices_list, function(mat) {
      t(scale(t(mat)))
    })
    
    rf <- JRF::JRF(
      X = norm_list,
      genes.name = rownames(norm_list[[1]]),
      ntree = 500,
      mtry = round(sqrt(nrow(norm_list[[1]]) - 1))
    )
    
    # Split JRF output into a list of data frames by 'importance' columns
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
  
  # For other methods: run in parallel for each matrix
  n_matrices <- length(count_matrices_list)
  nCores_outer <- min(total_cores, n_matrices)
  nCores_inner <- max(floor(total_cores / nCores_outer), 1)
  param_outer <- if (method == "GRNBoost2") BiocParallel::SerialParam() else BiocParallel::MulticoreParam(workers = nCores_outer, RNGseed = seed)
  
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
