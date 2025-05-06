#plotg
#' @keywords internal
#' @noRd
.create_igraph_plot <- function(mat, index) {
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = NULL, diag = FALSE)
  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
  
  if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) return(NULL)

  title <- paste("Graph", index, "\nNodes:", igraph::vcount(g), "Edges:", igraph::ecount(g))
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(color = "gray", width = 0.5) +
    ggraph::geom_node_point(color = "steelblue", size = 3) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )
}

#compare_consensus
#' @keywords internal
#' @noRd

.edge_to_str <- function(edge_list) {
  apply(edge_list, 1, function(e) paste(sort(e), collapse = "-"))
}
#' @keywords internal
#' @noRd

.plot_tp_fn_graph <- function(graph_ref, edge_colors, TP_label, FN_label) {
  graph_clean <- igraph::delete_vertices(graph_ref, igraph::V(graph_ref)[igraph::degree(graph_ref) == 0])

  ggraph(graph_clean, layout = "fr") +
    ggraph::geom_edge_link(aes(color = I(edge_colors)), width = 0.7) +
    ggraph::geom_node_point(color = "steelblue", size = 1.5) +
    ggplot2::labs(
      title = paste("Reference Graph\n", TP_label, ":", sum(edge_colors == "red"),
                    FN_label, ":", sum(edge_colors == "blue"))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}
#' @keywords internal
#' @noRd

.plot_fp_graph <- function(fp_edges_str, FP_label) {
  fp_pairs <- strsplit(fp_edges_str, "-")
  fp_mat <- do.call(rbind, Filter(function(x) length(x) == 2, fp_pairs))
  if (is.null(fp_mat) || nrow(fp_mat) == 0) return(NULL)

  graph_fp <- igraph::graph_from_edgelist(fp_mat, directed = FALSE)
  graph_fp <- igraph::delete_vertices(graph_fp, igraph::V(graph_fp)[igraph::degree(graph_fp) == 0])

  ggraph(graph_fp, layout = "fr") +
    ggraph::geom_edge_link(color = "purple", width = 1) +
    ggraph::geom_node_point(color = "steelblue", size = 2) +
    ggplot2::labs(title = paste(FP_label, ":", nrow(fp_mat))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#zinbsim
#' @keywords internal
#' @noRd
.compute_confusion_metrics <- function(pred_vec, gt_vec, index) {
  # — Argument checks —
  if (length(pred_vec) != length(gt_vec)) {
    stop("`pred_vec` and `gt_vec` must have the same length.")
  }
  if (!all(gt_vec %in% c(0, 1))) {
    stop("`gt_vec` must be binary (0/1).")
  }
  
  # force to integer 0/1
  pred_bin <- as.integer(pred_vec == 1)
  gt_bin   <- as.integer(gt_vec   == 1)
  
  TP <- sum(pred_bin & gt_bin)
  TN <- sum(!pred_bin & !gt_bin)
  FP <- sum(pred_bin & !gt_bin)
  FN <- sum(!pred_bin & gt_bin)
  
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
  Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  F1 <- ifelse((Precision + TPR) > 0, 2 * (Precision * TPR) / (Precision + TPR), 0)
  
  denominator <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                        as.numeric(TN + FP) * as.numeric(TN + FN))
  MCC <- ifelse(denominator > 0, (TP * TN - FP * FN) / denominator, 0)
  
  data.frame(
    Predicted_Matrix = paste("Matrix", index),
    TP, TN, FP, FN, TPR, FPR, Precision, F1, MCC
  )
}

#pscores
#' @keywords internal
#' @noRd

.plot_metrics_radar <- function(stats_df, metric_cols) {
  if (!all(metric_cols %in% colnames(stats_df))) {
    stop("Some metric columns are missing from stats_df.")
  }
  
  # pull raw and coerce
  score_data <- stats_df[, metric_cols, drop = FALSE]
  score_data <- as.data.frame(
    lapply(score_data, function(x) as.numeric(as.character(x)))
  )
  
  # clamp negative MCC → 0, with message
  if ("MCC" %in% colnames(score_data)) {
    neg <- which(score_data$MCC < 0)
    if (length(neg)) {
      message("Found ", length(neg),
              " negative MCC value(s); setting them to 0.")
      score_data$MCC[neg] <- 0
    }
  }
  
  # drop any non‐finite columns
  ok_cols <- colSums(is.finite(as.matrix(score_data))) > 0
  score_data <- score_data[, ok_cols, drop = FALSE]
  
  if (nrow(score_data) == 0 || ncol(score_data) < 2) {
    warning("Radar plot skipped: not enough valid metrics.")
    return(list(data = NULL, plot = NULL))
  }
  
  # decide axis max: if the highest value ≤ 0.5, use 0.5; else 1
  overall_max <- max(as.matrix(score_data), na.rm = TRUE)
  axis_max    <- if (overall_max <= 0.5) 0.5 else 1
  axis_min    <- 0
  
  # build the data.frame with max/min + observations
  max_row   <- rep(axis_max, ncol(score_data))
  min_row   <- rep(axis_min, ncol(score_data))
  plot_data <- rbind(max_row, min_row, score_data)
  
  rownames(plot_data) <- c("Max", "Min", stats_df$Predicted_Matrix)
  axis_labels <- pretty(c(axis_min, axis_max), n = 5)
  cols        <- grDevices::rainbow(nrow(score_data))
  
  graphics::par(mar = c(2, 2, 2, 2))
  fmsb::radarchart(
    data.frame(plot_data),
    axistype    = 2,
    pcol        = cols,
    plty        = 1,
    plwd        = 2,
    cglcol      = "grey",
    caxislabels = axis_labels,
    vlcex       = 1.1
  )
  graphics::legend(
    "topright",
    legend = stats_df$Predicted_Matrix,
    col    = cols,
    lty    = 1,
    lwd    = 2
  )
  
  return(list(data = plot_data, plot = grDevices::recordPlot()))
}
#' @keywords internal
#' @noRd
.plot_metrics_radar <- function(stats_df, metric_cols) {
  if (!all(metric_cols %in% colnames(stats_df))) {
    stop("Some metric columns are missing from stats_df.")
  }
  
  # pull raw and coerce
  score_data <- stats_df[, metric_cols, drop = FALSE]
  score_data <- as.data.frame(
    lapply(score_data, function(x) as.numeric(as.character(x)))
  )
  
  # clamp negative MCC → 0, with message
  if ("MCC" %in% colnames(score_data)) {
    neg <- which(score_data$MCC < 0)
    if (length(neg)) {
      message("Found ", length(neg),
              " negative MCC value(s); setting them to 0.")
      score_data$MCC[neg] <- 0
    }
  }
  
  # drop any non‐finite columns
  ok_cols <- colSums(is.finite(as.matrix(score_data))) > 0
  score_data <- score_data[, ok_cols, drop = FALSE]
  
  if (nrow(score_data) == 0 || ncol(score_data) < 2) {
    warning("Radar plot skipped: not enough valid metrics.")
    return(list(data = NULL, plot = NULL))
  }
  
  # decide axis max: if the highest value ≤ 0.5, use 0.5; else 1
  overall_max <- max(as.matrix(score_data), na.rm = TRUE)
  axis_max    <- if (overall_max <= 0.5) 0.5 else 1
  axis_min    <- 0
  
  # build the data.frame with max/min + observations
  max_row   <- rep(axis_max, ncol(score_data))
  min_row   <- rep(axis_min, ncol(score_data))
  plot_data <- rbind(max_row, min_row, score_data)
  
  rownames(plot_data) <- c("Max", "Min", stats_df$Predicted_Matrix)
  axis_labels <- pretty(c(axis_min, axis_max), n = 5)
  cols        <- grDevices::rainbow(nrow(score_data))
  
  graphics::par(mar = c(2, 2, 2, 2))
  fmsb::radarchart(
    data.frame(plot_data),
    axistype    = 2,
    pcol        = cols,
    plty        = 1,
    plwd        = 2,
    cglcol      = "grey",
    caxislabels = axis_labels,
    vlcex       = 1.1
  )
  graphics::legend(
    "topright",
    legend = stats_df$Predicted_Matrix,
    col    = cols,
    lty    = 1,
    lwd    = 2
  )
  
  return(list(data = plot_data, plot = grDevices::recordPlot()))
}

#earlyj
#' @keywords internal
#' @noRd

.merge_matrix_list <- function(input_list, rowg) {
  lapply(seq_along(input_list), function(i) {
    mat <- input_list[[i]]
    if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
      stop("Each matrix must be a non-empty numeric matrix.")
    }
    if (!rowg) mat <- t(mat)
    if (is.null(colnames(mat))) {
      colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
    }
    colnames(mat) <- paste0(colnames(mat), "-m", i)
    mat
  }) |> do.call(cbind, args = _)
}
#' @keywords internal
#' @noRd

.merge_seurat_list <- function(input_list) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required but not installed.")
  common_features <- Reduce(intersect, lapply(input_list, rownames))
  if (length(common_features) == 0) stop("No common features among Seurat objects.")

  modified <- lapply(seq_along(input_list), function(i) {
    obj <- Seurat::subset(input_list[[i]], features = common_features)
    Seurat::RenameCells(obj, new.names = paste0(Seurat::Cells(obj), "-m", i))
  })

  do.call(Seurat::merge, modified)
}
#' @keywords internal
#' @noRd

.merge_sce_list <- function(input_list) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package is required but not installed.")
  }

  common_features <- Reduce(intersect, lapply(input_list, rownames))
  if (length(common_features) == 0) stop("No common features among SCE objects.")

  modified <- lapply(seq_along(input_list), function(i) {
    obj <- input_list[[i]][common_features, ]
    colnames(obj) <- paste0(colnames(obj), "-m", i)
    obj
  })

  do.call(cbind, modified)
}

#cutoff_adjacency
#' @keywords internal
#' @noRd

.shuffle_matrix_rows <- function(mat) {
  shuffled <- t(apply(mat, 1, sample))
  rownames(shuffled) <- rownames(mat)
  colnames(shuffled) <- colnames(mat)  
  return(shuffled)
}

.run_network_on_shuffled <- function(mat, method, grnboost_modules, weight_function, quantile_threshold) {
  shuffled <- .shuffle_matrix_rows(mat)
  inferred <- infer_networks(list(shuffled), method = method, grnboost_modules = grnboost_modules)
  adjm <- generate_adjacency(inferred)
  symm <- symmetrize(adjm, weight_function = weight_function)[[1]]
  quantile(symm[upper.tri(symm)], quantile_threshold, names = FALSE)
}
#' @keywords internal
#' @noRd

.aggregate_cutoffs <- function(results, n_matrices) {
  percentile_values_by_matrix <- vector("list", n_matrices)
  for (res in results) {
    mat_idx <- res$matrix_idx
    percentile_values_by_matrix[[mat_idx]] <- c(percentile_values_by_matrix[[mat_idx]], res$q_value)
  }
  percentile_values_by_matrix
}
#' @keywords internal
#' @noRd

.binarize_adjacency <- function(weighted_list, cutoffs, method, debug) {
  lapply(seq_along(weighted_list), function(i) {
    avg_cutoff <- mean(cutoffs[[i]])
    if (debug) message(sprintf("[Method: %s] Matrix %d → Cutoff = %.5f", method, i, avg_cutoff))
    ifelse(weighted_list[[i]] > avg_cutoff, 1, 0)
  })
}

#plotROC
#' @keywords internal
#' @noRd

.prepare_prediction_vectors <- function(pred_matrix, ground_truth) {
  pred_matrix <- as.data.frame(pred_matrix)
  pred_matrix <- pred_matrix[rownames(ground_truth), colnames(ground_truth)]
  as.vector(as.matrix(pred_matrix)[upper.tri(pred_matrix)])
}
#' @keywords internal
#' @noRd

.compute_binary_roc_point <- function(pred_vec, truth_vec) {
  tp <- sum(pred_vec == 1 & truth_vec == 1)
  fp <- sum(pred_vec == 1 & truth_vec == 0)
  fn <- sum(pred_vec == 0 & truth_vec == 1)
  tn <- sum(pred_vec == 0 & truth_vec == 0)

  TPR <- tp / (tp + fn)
  FPR <- fp / (fp + tn)
  list(FPR = FPR, TPR = TPR)
}
#' @keywords internal
#' @noRd

.compute_weighted_roc_curve <- function(pred_vec, truth_vec, label) {
  roc_obj <- pROC::roc(truth_vec, pred_vec, quiet = TRUE)
  auc <- round(roc_obj$auc, 2)

  df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Matrix = paste0(label, " (AUC=", auc, ")")
  )
  list(df = df, auc = auc)
}
#' @keywords internal
#' @noRd

.plot_roc_curve <- function(roc_data, binary_points = NULL, title = "") {
  total_matrices <- length(unique(roc_data$Matrix))
  colors <- scales::hue_pal()(total_matrices)

  p <- ggplot() +
    labs(title = title, x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    geom_line(data = roc_data, aes(x = FPR, y = TPR, color = Matrix), size = 1.2)

  if (!is.null(binary_points)) {
    p <- p + geom_point(data = binary_points, aes(x = FPR, y = TPR), size = 2, color = "blue")
  }

  p + scale_color_manual(values = colors)
}


#selgene
#' @keywords internal
#' @noRd

.extract_expression <- function(object, assay = NULL) {
  if (inherits(object, "Seurat")) {
    assay_name <- Seurat::DefaultAssay(object)
    seurat_assay <- object[[assay_name]]
    slots_avail <- methods::slotNames(seurat_assay)

    if (!"data" %in% slots_avail) {
      stop("Assay '", assay_name, "' has no 'data' slot. Available slots: ",
           paste(slots_avail, collapse = ", "))
    }
    message("Using Seurat assay '", assay_name, "' slot 'data' (log-normalized).")
    return(seurat_assay@data)
  }

  if (inherits(object, "SingleCellExperiment")) {
    available_assays <- SummarizedExperiment::assayNames(object)
    assay_to_use <- if (!is.null(assay)) assay else "logcounts"

    if (!assay_to_use %in% available_assays) {
      stop("Requested assay '", assay_to_use, "' not found. Available assays: ",
           paste(available_assays, collapse = ", "))
    }

    message("Using SCE assay '", assay_to_use, "' (assumed log-normalized).")
    return(SummarizedExperiment::assay(object, assay_to_use))
  }

  if (is.matrix(object)) {
    return(object)
  }

  stop("Input must be a Seurat, SingleCellExperiment, or matrix.")
}
#' @keywords internal
#' @noRd

.filter_by_cell_type <- function(expr, object, cell_type, cell_type_col) {
  if (inherits(object, "Seurat")) {
    meta <- object@meta.data
  } else if (inherits(object, "SingleCellExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(object))
  } else {
    stop("'cell_type' filtering is not supported for raw matrices.")
  }

  if (!cell_type_col %in% colnames(meta)) {
    stop("Metadata must contain column '", cell_type_col, "'.")
  }

  keep_cells <- rownames(meta)[meta[[cell_type_col]] == cell_type]
  message("Subsetted to ", length(keep_cells), " cells where ", cell_type_col, " = '", cell_type, "'.")
  return(expr[, colnames(expr) %in% keep_cells, drop = FALSE])
}
#' @keywords internal
#' @noRd

.filter_genes <- function(expr, remove_mt, remove_rib) {
  gene_names <- rownames(expr)
  keep_genes <- rep(TRUE, length(gene_names))

  if (remove_mt) {
    keep_genes <- keep_genes & !grepl("^MT-", gene_names, ignore.case = TRUE)
    message("Removed mitochondrial genes matching '^MT-'.")
  }

  if (remove_rib) {
    keep_genes <- keep_genes & !grepl("^RP[SL]", gene_names, ignore.case = TRUE)
    message("Removed ribosomal genes matching '^RP[SL]'.")
  }

  expr[keep_genes, , drop = FALSE]
}
#' @keywords internal
#' @noRd

.select_top_genes <- function(expr, top_n) {
  expr <- as.matrix(expr)
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]

  avg_expression <- rowMeans(expr, na.rm = TRUE)
  sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
  head(sorted_genes, top_n)
}


#Infer_ networks
#' @keywords internal
#' @noRd

.convert_counts_list <- function(count_matrices_list) {
  lapply(count_matrices_list, function(obj) {
    if (inherits(obj, "Seurat")) {
      as.matrix(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
    } else if (inherits(obj, "SingleCellExperiment")) {
      as.matrix(SummarizedExperiment::assay(obj, "counts"))
    } else {
      as.matrix(obj)
    }
  })
}
#' @keywords internal
#' @noRd

.run_genie3 <- function(mat, nCores) {
  adj <- GENIE3::GENIE3(mat, nCores = nCores)
  GENIE3::getLinkList(adj)
}
#' @keywords internal
#' @noRd

.run_zilgm <- function(mat, adjm, nCores) {
  lambda_max <- ZILGM::find_lammax(t(mat))
  lambda_seq <- exp(seq(log(lambda_max), log(1e-4 * lambda_max), length.out = 50))
  fit <- ZILGM::zilgm(
    X = t(mat), lambda = lambda_seq, family = "NBII",
    update_type = "IRLS", do_boot = TRUE, boot_num = 10,
    sym = "OR", nCores = nCores
  )
  adj <- fit$network[[fit$opt_index]]
  dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
  adj
}
#' @keywords internal
#' @noRd

.run_jrf <- function(norm_list, nCores) {
  clust <- parallel::makeCluster(nCores)
  on.exit(parallel::stopCluster(clust), add = TRUE)
  doParallel::registerDoParallel(clust)

  rf <- JRF::JRF(
    X = norm_list,
    genes.name = rownames(norm_list[[1]]),
    ntree = 500,
    mtry = round(sqrt(nrow(norm_list[[1]]) - 1))
  )

  importance_columns <- grep("importance", names(rf), value = TRUE)
  lapply(importance_columns, function(col) {
    df <- rf[, c("gene1", "gene2", col)]
    names(df)[3] <- col
    df
  })
}
#' @keywords internal
#' @noRd

.run_parallel_networks <- function(count_matrices_list, method, nCores, adjm, grnboost_modules) {
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
      return(result_r)
    } else if (method == "PCzinb") {
      adj <- learn2count::PCzinb(t(mat), method = "zinb1", maxcard = 2)
      dimnames(adj) <- if (is.null(adjm)) list(rownames(mat), rownames(mat)) else dimnames(adjm)
      return(adj)
    }
  }, BPPARAM = param_outer)
}


#This is community similarity
#' @keywords internal
#' @noRd

.compute_topo_metrics <- function(graph, comm) {
  c(
    Modularity = igraph::modularity(graph, comm),
    Communities = length(unique(comm)),
    Density = igraph::edge_density(graph),
    Transitivity = igraph::transitivity(graph)
  )
}
#' @keywords internal
#' @noRd

.compare_communities <- function(control_comm, pred_comm) {
  c(
    VI = igraph::compare(control_comm, pred_comm, method = "vi"),
    NMI = igraph::compare(control_comm, pred_comm, method = "nmi"),
    ARI = igraph::compare(control_comm, pred_comm, method = "adjusted.rand")
  )
}
#' @keywords internal
#' @noRd

.plot_radar_communities <- function(comm_df) {
  max_val <- ceiling(max(comm_df, na.rm = TRUE))
  axis_steps <- pretty(c(0, max_val), n = 5)
  radar_comm <- rbind(rep(max_val, ncol(comm_df)), rep(0, ncol(comm_df)), comm_df)

  colors <- grDevices::rainbow(nrow(comm_df))
  graphics::par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

  fmsb::radarchart(radar_comm,
    axistype = 2, pcol = colors, plwd = 2, plty = 1,
    cglcol = "grey", axislabcol = "black",
    caxislabels = axis_steps, vlcex = 1.1,
    title = "Community Similarity Metrics"
  )

  graphics::legend("topright", legend = rownames(comm_df), col = colors, lty = 1, lwd = 2)
}
#' @keywords internal
#' @noRd

.plot_topo_barplots <- function(topo_df, control_topo) {
  for (i in seq_len(nrow(topo_df))) {
    pred_name <- rownames(topo_df)[i]
    pred_topo <- topo_df[i, ]

    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))
    for (metric in colnames(pred_topo)) {
      graphics::barplot(
        height = c(control_topo[[metric]], pred_topo[[metric]]),
        names.arg = c("Control", "Predicted"),
        main = paste0(metric, " Comparison\n", pred_name),
        ylab = metric,
        col = c("lightblue", "salmon")
      )
    }
  }
}

#This is edge_mining.R
#' @keywords internal
#' @noRd

.identify_edges <- function(predicted, ground_truth, query_edge_types) {
  indices <- which(((predicted == 1) | (ground_truth == 1)) & upper.tri(predicted), arr.ind = TRUE)

  if (nrow(indices) == 0) return(NULL)

  gene_pairs <- data.frame(
    gene1 = rownames(predicted)[indices[, "row"]],
    gene2 = colnames(predicted)[indices[, "col"]],
    stringsAsFactors = FALSE
  )

  gene_pairs$edge_type <- ifelse(predicted[indices] == 1 & ground_truth[indices] == 1, "TP",
    ifelse(predicted[indices] == 1 & ground_truth[indices] == 0, "FP", "FN")
  )

  gene_pairs <- gene_pairs[gene_pairs$edge_type %in% query_edge_types, , drop = FALSE]

  if (nrow(gene_pairs) == 0) return(NULL)
  return(gene_pairs)
}
#' @keywords internal
#' @noRd

.safe_query_pubmed <- function(gene1, gene2, query_field, delay, max_retries) {
  query <- paste0(gene1, "[", query_field, "] AND ", gene2, "[", query_field, "]")

  for (attempt in seq_len(max_retries)) {
    result <- tryCatch(
      {
        search_res <- rentrez::entrez_search(db = "pubmed", term = query, retmax = 100)
        Sys.sleep(delay)
        list(
          pubmed_hits = as.numeric(search_res$count),
          PMIDs = if (length(search_res$ids) > 0) paste(search_res$ids, collapse = ",") else NA_character_
        )
      },
      error = function(e) NULL
    )
    if (!is.null(result)) return(result)
    Sys.sleep(delay)
  }

  return(list(pubmed_hits = NA_integer_, PMIDs = NA_character_))
}
#' @keywords internal
#' @noRd

.query_edge_pairs <- function(gene_pairs, query_field, delay, max_retries, BPPARAM) {
  pubmed_info <- BiocParallel::bplapply(seq_len(nrow(gene_pairs)), function(j) {
    res <- .safe_query_pubmed(gene_pairs$gene1[j], gene_pairs$gene2[j], query_field, delay, max_retries)
    data.frame(pubmed_hits = res$pubmed_hits, PMIDs = res$PMIDs, stringsAsFactors = FALSE)
  }, BPPARAM = BPPARAM)

  pubmed_info <- do.call(rbind, pubmed_info)

  gene_pairs$pubmed_hits <- pubmed_info$pubmed_hits
  gene_pairs$PMIDs <- pubmed_info$PMIDs
  gene_pairs$query_status <- ifelse(is.na(gene_pairs$pubmed_hits), "error",
                             ifelse(gene_pairs$pubmed_hits == 0, "no_hits", "hits_found"))

  return(gene_pairs)
}

#This is stringdb
#' @keywords internal
#' @noRd

.map_genes_to_string <- function(string_db, genes) {
  df <- data.frame(genes, stringsAsFactors = FALSE)
  mapped <- string_db$map(df, "genes", removeUnmappedRows = FALSE)
  mapped_clean <- mapped[!is.na(mapped$STRING_id), ]
  unmapped <- setdiff(genes, mapped_clean$genes)
  return(list(mapped = mapped_clean, unmapped = unmapped))
}
#' @keywords internal
#' @noRd

.query_string_api <- function(mapped_ids, species, required_score) {
  base_url <- "https://string-db.org/api/json/network"
  identifiers_str <- paste(mapped_ids, collapse = "\n")

  res <- httr::POST(
    url = base_url,
    body = list(
      identifiers    = identifiers_str,
      species        = species,
      required_score = required_score,
      network_type   = "physical"
    ),
    encode = "form"
  )

  if (res$status_code != 200) {
    stop("STRING API query failed. Status code: ", res$status_code)
  }

  interactions <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))

  if (!is.data.frame(interactions) || nrow(interactions) == 0) {
    return(NULL)
  }

  return(interactions)
}
#' @keywords internal
#' @noRd

.build_adjacency_matrices <- function(interactions, mapped_genes, genes, keep_all_genes) {
  interactions$interaction_score <- interactions$score

  id_to_gene <- setNames(mapped_genes$genes, mapped_genes$STRING_id)
  interactions$gene_A <- id_to_gene[interactions$stringId_A]
  interactions$gene_B <- id_to_gene[interactions$stringId_B]

  interactions <- interactions[!is.na(interactions$gene_A) & !is.na(interactions$gene_B), ]

  final_gene_list <- if (keep_all_genes) genes else unique(c(interactions$gene_A, interactions$gene_B))
  p <- length(final_gene_list)

  weighted_mat <- matrix(0, nrow = p, ncol = p, dimnames = list(final_gene_list, final_gene_list))
  for (i in seq_len(nrow(interactions))) {
    a <- interactions$gene_A[i]
    b <- interactions$gene_B[i]
    s <- interactions$interaction_score[i]
    if (!is.na(a) && !is.na(b) && a %in% final_gene_list && b %in% final_gene_list) {
      weighted_mat[a, b] <- s
      weighted_mat[b, a] <- s
    }
  }

  binary_mat <- ifelse(weighted_mat > 0, 1, 0)

  return(list(weighted = weighted_mat, binary = binary_mat))
}
#' @keywords internal
#' @noRd

.zero_matrix_result <- function(genes) {
  zero_mat <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  list(weighted = zero_mat, binary = zero_mat)
}


# This is community_path
#' @keywords internal
#' @noRd

.detect_communities <- function(graph, methods) {
  if (length(methods) == 1) {
    best_method <- methods[1]
    best_communities <- robin::membershipCommunities(graph, method = best_method)
  } else if (length(methods) == 2) {
    res <- tryCatch(
      robin::robinCompare(graph, method1 = methods[1], method2 = methods[2]),
      error = function(e) stop("robinCompare failed: ", conditionMessage(e))
    )

    auc <- robin::robinAUC(res)
    if (!is.numeric(auc) || length(auc) != 2) stop("Unexpected AUC result from robinCompare.")

    if (auc[1] < auc[2]) {
      best_method <- methods[1]
      best_communities <- res$Communities1
    } else {
      best_method <- methods[2]
      best_communities <- res$Communities2
    }
  } else {
    stop("methods must be a character vector of length 1 or 2.")
  }

  return(list(best_method = best_method, best_communities = best_communities))
}
#' @keywords internal
#' @noRd

.plot_communities <- function(graph, best_method) {
  non_isolated_nodes <- igraph::degree(graph) > 0
  plot_graph <- igraph::induced_subgraph(graph, vids = igraph::V(graph)[non_isolated_nodes])
  num_communities <- length(unique(igraph::V(plot_graph)$community))

  colors <- if (num_communities <= 12) {
    RColorBrewer::brewer.pal(num_communities, "Set3")
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_communities)
  }

  plot_title <- paste0(
    "Community Structure (", best_method, ")\nNodes: ",
    igraph::vcount(plot_graph), " Edges: ", igraph::ecount(plot_graph)
  )

  g <- ggraph::ggraph(plot_graph, layout = "fr") +
    ggraph::geom_edge_link(color = "gray", width = 0.5) +
    ggraph::geom_node_point(ggplot2::aes(color = community), size = 3) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )

  plot(g)
}
#' @keywords internal
#' @noRd

.enrich_communities <- function(graph, non_isolated_nodes, pathway_db, genes_path) {
  pathway_results <- list()
  non_isolated_genes <- igraph::V(graph)$name[non_isolated_nodes]

  for (comm in unique(igraph::V(graph)$community)) {
    genes <- intersect(
      igraph::V(graph)$name[igraph::V(graph)$community == comm],
      non_isolated_genes
    )

    if (length(genes) < genes_path) {
      pathway_results[[as.character(comm)]] <- NULL
      next
    }

    entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
      keys = genes,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    entrez <- na.omit(entrez)

    if (length(entrez) >= genes_path) {
      enrich <- tryCatch(
        {
          switch(pathway_db,
            "KEGG" = clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa", keyType = "kegg"),
            "Reactome" = ReactomePA::enrichPathway(gene = entrez, organism = "human"),
            NULL
          )
        },
        error = function(e) {
          warning("Enrichment failed for community ", comm, ": ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(enrich) && nrow(enrich@result) > 0) {
        pathway_results[[as.character(comm)]] <- enrich
      } else {
        pathway_results[[as.character(comm)]] <- NULL
      }
    } else {
      pathway_results[[as.character(comm)]] <- NULL
    }
  }

  return(pathway_results)
}

